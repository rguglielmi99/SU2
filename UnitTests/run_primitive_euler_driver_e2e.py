#!/usr/bin/env python3
"""End-to-end SU2_CFD tests, using full-precision native restart output.

VTU is intentionally checked separately: SU2's native XML writer uses Float32.
No direct solver calls or geometry construction are used by the driver runs.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import struct
import subprocess
import time
import xml.etree.ElementTree as ET
import numpy as np
from primitive_euler_driver_cases import mesh, config, ROOT


def binary(path):
    data=path.read_bytes()
    magic,nvar,npoint,_,_=struct.unpack_from('=5i',data)
    if magic!=535532:raise ValueError('Not a native SU2 binary restart')
    fields=[data[20+33*k:20+33*(k+1)].split(b'\0')[0].decode() for k in range(nvar)]
    values=np.frombuffer(data,dtype=np.float64,count=nvar*npoint,offset=20+33*nvar).reshape(npoint,nvar)
    if not np.all(np.isfinite(values)):raise ValueError('Nonfinite native output')
    return dict(zip(fields,values.T))


def reference(family,n,index):
    piece=ET.parse(ROOT/f'Docs/primitive_euler_phase7_vortex/{family}_{n}_{index}.vtu').getroot().find('UnstructuredGrid/Piece')
    point={a.attrib['Name']:np.fromstring(a.text,sep=' ') for a in piece.findall('PointData/DataArray')}
    xy=np.fromstring(piece.find('Points/DataArray').text,sep=' ').reshape(-1,3)[:,:2]
    return xy,point


def write(path,rows):
    with path.open('w') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def reference_checks(exe,out):
    results=[]
    for family in ('regular','jitter'):
        n=16
        path=out/f'reference_{family}_{n}';path.mkdir(exist_ok=True)
        with (path/'run.log').open('w') as log:
            subprocess.run([str(exe.resolve()),family,str(n),'8',str(path/'reference.csv')],check=True,stdout=log,stderr=subprocess.STDOUT,
                           env={**os.environ,'OMP_NUM_THREADS':'1'})
        def geometry(file):
            values={}
            with file.open() as f:
                for row in csv.DictReader(f):
                    kind=row['kind'];i,j=int(row['first']),int(row['second'])
                    v=np.array([float(row[k]) for k in ('x','y','nx','ny','volume')])
                    if kind=='edge' and i>j:i,j=j,i;v[2:4]*=-1
                    key=(kind,i,j,round(v[0],13),round(v[1],13)) if kind=='boundary' else (kind,i,j)
                    if key in values:raise ValueError('Duplicate metric')
                    values[key]=v
            return values
        native=geometry(out/f'vortex_{family}_{n}/geometry.csv')
        api=geometry(path/'reference.csv.geometry.csv')
        if native.keys()!=api.keys():raise ValueError('Geometry identities differ')
        maxgeom=max(np.max(abs(native[k]-api[k])) for k in native)
        if maxgeom>1e-13:raise ValueError(f'Geometry divergence {maxgeom}')
        maxstate=0
        for index,t in enumerate((0.,.25,.5)):
            field=np.loadtxt(str(path/'reference.csv')+f'.snapshot_{t:.6f}')
            _,old=reference(family,n,index)
            original=np.column_stack([old[name+'_primitive'] for name in ('rho','u','v','p')])
            maxstate=max(maxstate,float(np.max(abs(field[:,4:8]-original))))
        if maxstate>2e-13:raise ValueError('Rebuilt phase-7 harness differs from archived results')
        results.append(dict(family=family,N=n,geometry_max_difference=maxgeom,rebuilt_harness_max_difference=maxstate,
                            metrics=len(native),boundary_count=sum(k[0]=='boundary' for k in native)))
    write(out/'geometry_reference.csv',results)


def run(exe,path,text,timeout=360):
    path.mkdir(parents=True,exist_ok=True)
    (path/'case.cfg').write_text(text)
    start=time.monotonic()
    with (path/'run.log').open('w') as log:
        result=subprocess.run([str(exe),'case.cfg'],cwd=path,stdout=log,stderr=subprocess.STDOUT,
            env={**os.environ,'OMP_NUM_THREADS':'1'},timeout=timeout)
    return result.returncode,time.monotonic()-start


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--pilot',action='store_true')
    p.add_argument('--negative-only',action='store_true')
    p.add_argument('--reference-only',action='store_true')
    p.add_argument('--reference-executable',type=Path)
    a=p.parse_args()
    exe=a.executable.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    if a.reference_only:
        if a.reference_executable is None:p.error('--reference-executable required')
        reference_checks(a.reference_executable,out)
        return
    comparisons,errors,uniform,execution=[] ,[],[],[]
    if not a.negative_only:
        for family in ('regular','jitter'):
            for kind in (['uniform'] if a.pilot else ['uniform','vortex']):
                for n in ([8] if kind=='uniform' or a.pilot else [8,16,32,64]):
                    path=out/f'{kind}_{family}_{n}'
                    (path).mkdir(exist_ok=True)
                    (path/'mesh.su2').write_text(mesh(family,n))
                    code,seconds=run(exe,path,config('mesh.su2',n,kind=='uniform'))
                    if code:raise RuntimeError(f'{path}: failed ({code}); inspect run.log')
                    histories=list(path.glob('history*.csv'))
                    if len(histories)!=1:raise ValueError('Missing/ambiguous native history')
                    with histories[0].open() as stream:
                        reader=csv.DictReader(stream,skipinitialspace=True)
                        reader.fieldnames=[name.strip() for name in reader.fieldnames]
                        history=list(reader)
                    time_at_step={int(float(r['Time_Iter'])):float(r['Cur_Time']) for r in history}
                    with (path/'geometry.csv').open() as stream:
                        volumes={int(r['first']):float(r['volume']) for r in csv.DictReader(stream) if r['kind']=='node'}
                    original_xy,_=reference(family,n,0)
                    point_ids={tuple(x):i for i,x in enumerate(original_xy)}
                    records=[]
                    for f in sorted(path.glob('primitive_state*.dat')):
                        b=binary(f)
                        step_index=int(f.stem.rsplit('_',1)[1])
                        t=time_at_step[step_index]
                        # Native restart keeps only coordinates+SOLUTION, in Float64.
                        # Associate the checkpoint with its native history index and
                        # match coordinates to global IDs; do not infer IDs from RCM.
                        b['Global_ID']=np.array([point_ids[x] for x in zip(b['x'],b['y'])])
                        b['Dual_Volume']=np.array([volumes[i] for i in b['Global_ID']])
                        records.append((t,b,f))
                    if len(records)!=3 or len({r[0] for r in records})!=3:
                        raise ValueError(f'Wrong snapshots in {path}: {[r[0] for r in records]}')
                    final=.08 if kind=='uniform' else .5
                    np.testing.assert_allclose(sorted(r[0] for r in records),[0,final/2,final],rtol=0,atol=64*np.finfo(float).eps)
                    expected_steps=8 if kind=='uniform' else 2*n*n
                    if len(history)!=expected_steps+1:raise ValueError('History does not contain initial state plus all accepted steps')
                    np.testing.assert_array_equal([int(float(r['Time_Iter'])) for r in history],np.arange(expected_steps+1))
                    np.testing.assert_allclose(float(history[-1]['Cur_Time']),final,rtol=0,atol=64*np.finfo(float).eps)
                    execution.append(dict(kind=kind,family=family,N=n,seconds=seconds,snapshots=len(records),returncode=code))
                    for index,(t,b,f) in enumerate(sorted(records,key=lambda r:r[0])):
                        ids=b['Global_ID'].astype(int)
                        order=np.argsort(ids)
                        np.testing.assert_array_equal(ids[order],np.arange((n+1)**2))
                        w=np.column_stack([b[k][order] for k in ('Density','Velocity_x','Velocity_y','Pressure')])
                        volume=b['Dual_Volume'][order]
                        xy,ref=reference(family,n,index if kind=='vortex' else 0)
                        np.testing.assert_array_equal(np.column_stack((b['x'][order],b['y'][order])),xy)
                        np.testing.assert_allclose(volume,ref['dual_volume'],rtol=2e-14,atol=2e-15)
                        for k,name in enumerate(('rho','u','v','p')):
                            if kind=='uniform':
                                e=w[:,k]-(2,.5,.3,3)[k]
                                uniform.append(dict(family=family,N=n,time=t,component=k,max_error=float(np.max(abs(e)))))
                                if np.max(abs(e))>2e-13:raise ValueError('Uniform state changed')
                            else:
                                e=w[:,k]-ref[name+'_primitive']
                                row=dict(family=family,N=n,time=t,component=k,l1=np.sum(volume*abs(e))/volume.sum(),
                                    l2=np.sqrt(np.sum(volume*e*e)/volume.sum()),linf=np.max(abs(e)))
                                comparisons.append(row)
                                if row['linf']>2e-13:raise ValueError(f'Driver/harness diverged: {row}')
                                ea=w[:,k]-ref[name+'_exact']
                                errors.append(dict(family=family,N=n,time=t,component=k,l1=np.sum(volume*abs(ea))/volume.sum(),
                                    l2=np.sqrt(np.sum(volume*ea*ea)/volume.sum()),linf=np.max(abs(ea))))
                    print(kind,family,n,'passed',seconds,flush=True)
        write(out/'execution.csv',execution)
        write(out/'uniform.csv',uniform)
        if comparisons:
            write(out/'equivalence.csv',comparisons)
            previous={}
            for row in sorted(errors,key=lambda r:r['N']):
                key=row['family'],row['time'],row['component']
                old=previous.get(key)
                row['p_l2']=np.log2(old['l2']/row['l2']) if old and min(old['l2'],row['l2'])>1e-12 else ''
                previous[key]=row
            write(out/'errors.csv',errors)
    if a.pilot:return
    base=config('mesh.su2',8,True)
    variants={'implicit':('TIME_DISCRE_FLOW=EULER_EXPLICIT','TIME_DISCRE_FLOW=EULER_IMPLICIT','EULER_EXPLICIT'),
        'multigrid':('MGLEVEL=0','MGLEVEL=2','multigrid'),
        'dual_time':('TIME_MARCHING=TIME_STEPPING','TIME_MARCHING=DUAL_TIME_STEPPING-1ST_ORDER','TIME_STEPPING'),
        'unknown_marker':('MARKER_CUSTOM=(boundary)','MARKER_CUSTOM=(wrong)','marker'),
        'unused_marker':('MARKER_CUSTOM=(boundary)','MARKER_CUSTOM=(boundary,unused)','marker'),
        'wall':('MARKER_CUSTOM=(boundary)','MARKER_EULER=(boundary)','supports only'),
        'symmetry':('MARKER_CUSTOM=(boundary)','MARKER_SYM=(boundary)','supports only'),
        'bad_mode':('PRIMITIVE_INITIALIZATION=UNIFORM','PRIMITIVE_INITIALIZATION=QUADRATIC','UNIFORM'),
        'bad_representation':('PRIMITIVE_REPRESENTATION=NODE','PRIMITIVE_REPRESENTATION=DUAL_MEAN','NODE'),
        'bad_final_time':('MAX_TIME=0.080000000000000002','MAX_TIME=0.085','integer multiple'),
        'motion':('MGLEVEL=0','MGLEVEL=0\nGRID_MOVEMENT=RIGID_MOTION','static'),
        'cfl_conflict':('UNST_CFL_NUMBER=0','UNST_CFL_NUMBER=0.1','TIME_STEP=0')}
    variants['zero_output_frequency']=('OUTPUT_WRT_FREQ=4','OUTPUT_WRT_FREQ=0','must be positive')
    variants['native_mms']=('MGLEVEL=0','MGLEVEL=0\nKIND_VERIFICATION_SOLUTION=MMS_NS_UNIT_QUAD','KIND_VERIFICATION_SOLUTION')
    negative=[]
    for name,(find,replace,expected) in variants.items():
        path=out/('negative_'+name);path.mkdir(exist_ok=True)
        (path/'mesh.su2').write_text(mesh('regular',8))
        if find not in base:raise ValueError(find)
        code,seconds=run(exe,path,base.replace(find,replace))
        log=(path/'run.log').read_text()
        if code==0 or expected.lower() not in log.lower():raise ValueError(f'Negative case {name} failed: {log[-3000:]}')
        negative.append(dict(case=name,returncode=code,expected=expected,seconds=seconds))
    for name,mesh_text,expected in [('3d',mesh('regular',8).replace('NDIME= 2','NDIME= 3'),'2D'),
            ('quadrilateral','NDIME=2\nNELEM=1\n9 0 1 2 3\nNPOIN=4\n0 0 0\n1 0 1\n1 1 2\n0 1 3\nNMARK=1\nMARKER_TAG=boundary\nMARKER_ELEMS=4\n3 0 1\n3 1 2\n3 2 3\n3 3 0\n','triang')]:
        path=out/('negative_'+name);path.mkdir(exist_ok=True);(path/'mesh.su2').write_text(mesh_text)
        code,seconds=run(exe,path,base)
        if code==0 or expected.lower() not in (path/'run.log').read_text().lower():raise ValueError(name)
        negative.append(dict(case=name,returncode=code,expected=expected,seconds=seconds))
    write(out/'negative.csv',negative)


if __name__=='__main__':main()
