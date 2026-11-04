#!/usr/bin/env python3
"""Collect BC experiments, observed orders, temporal checks and reflection fits."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import numpy as np

PATTERN=re.compile(r'^(regular|jitter)_(.+)_(\d+)_(\d+)\.csv$')


def read(path):
    with path.open() as stream:return list(csv.DictReader(stream))


def write(path,rows):
    if not rows:raise ValueError('Empty result: '+str(path))
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def samples(directory,acoustic):
    result=[]
    for path in sorted(directory.glob('*.csv')):
        match=PATTERN.match(path.name)
        if not match or match[2].startswith('acoustic_')!=acoustic:continue
        rows=read(path)
        final=6 if acoustic else .5
        if not rows or max(float(r['time']) for r in rows)!=final:
            raise ValueError('Incomplete experiment: '+str(path))
        for r in rows:
            if any(not math.isfinite(float(r[k])) for k in ('l1','l2','linf','min_rho','min_p')):
                raise ValueError('Nonfinite norm')
            if min(float(r['min_rho']),float(r['min_p']))<=0:raise ValueError('Nonpositive state')
        result.extend(rows)
    return result


def global_errors(rows):
    final=[dict(r) for r in rows if float(r['time'])>0 and r['region']=='0']
    previous={}
    for r in sorted(final,key=lambda row:int(row['N'])):
        key=tuple(r[k] for k in ('family','mode','divisor','kind','component','amplitude'))
        old=previous.get(key);e=float(r['l2'])
        r['order_l2']=''
        if old and min(float(old['l2']),e)>1e-13:
            r['order_l2']=math.log(float(old['l2'])/e)/math.log(int(r['N'])/int(old['N']))
        previous[key]=r
    return final


def regional(rows):
    groups={}
    for r in rows:
        if r['kind'] not in ('P','tau_P'):continue
        if r['kind']=='P' and float(r['time'])==0:continue
        key=tuple(r[k] for k in ('family','N','mode','divisor','time','kind','region'))
        groups.setdefault(key,[]).append(r)
    result=[]
    for key,values in groups.items():
        result.append(dict(zip(('family','N','mode','divisor','time','kind','region'),key),
            l1_vector_one=sum(float(r['l1']) for r in values),
            l2_vector_euclidean=math.sqrt(sum(float(r['l2'])**2 for r in values)),
            linf_vector_max=max(float(r['linf']) for r in values)))
    return result


def reflections(directory):
    result=[]
    for path in sorted(directory.glob('*.csv.waves.csv')):
        match=PATTERN.match(path.name.removesuffix('.waves.csv'))
        if not match:continue
        normpath=path.with_name(path.name.removesuffix('.waves.csv'))
        norms=read(normpath);amplitude=float(norms[-1]['amplitude']);t=float(norms[-1]['time'])
        data=read(path);mode=match[2];u=0 if mode=='acoustic_wall' else .4;c=math.sqrt(1.4)
        left=mode in ('acoustic_total','acoustic_mass');direction=-1 if left else 1
        boundary=-5 if left else 5;center=-2.5 if left else 2.5
        x=np.array([float(r['x']) for r in data]);v=np.array([float(r['volume']) for r in data])
        delay=(x-boundary)/(u-direction*c)
        arg=boundary-center-(u+direction*c)*(t-delay)
        shape=np.where((abs(arg)<1.5)&(delay>=0)&(delay<=t),(1-arg*arg/2.25)**4,0)
        denominator=np.sum(v*shape*shape)
        if denominator<=0:raise ValueError('Unresolved reflection window')
        expected={'acoustic_far':0,'acoustic_wall':1,'acoustic_mass':1,'acoustic_outlet':-1,
                  'acoustic_total':-(1-u/c)/(1+u/c)}[mode]
        for solver in ('P','C'):
            pressure=np.array([float(r['p_'+solver])-1 for r in data])
            fit=float(np.sum(v*pressure*shape)/(amplitude*denominator))
            energy=float(np.sqrt(np.sum(v*pressure*pressure)/(amplitude*amplitude*denominator)))
            result.append(dict(family=match[1],mode=mode,N=int(match[3]),divisor=int(match[4]),
                solver=solver,amplitude=amplitude,R_linear=expected,R_fit=fit,
                R_fit_error=abs(fit-expected),pressure_energy_ratio=energy))
    return result


def temporal(rows):
    errors=global_errors(rows);index={tuple(r[k] for k in ('family','N','mode','kind','component','divisor')):r for r in errors}
    result=[]
    for r in errors:
        divisor=int(r['divisor']);key=tuple(r[k] for k in ('family','N','mode','kind','component'))
        fine=index.get(key+(str(2*divisor),))
        if not fine:continue
        e1=float(r['l2']);e2=float(fine['l2']);relative=abs(e1-e2)/max(e1,e2) if max(e1,e2)>1e-12 else 0
        result.append(dict(zip(('family','N','mode','kind','component'),key),divisor=divisor,
            error_dt=e1,error_half_dt=e2,relative_norm_change=relative,below_roundoff_floor=int(max(e1,e2)<=1e-12),
            passes_five_percent=int(relative<=.05)))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('contact','acoustic','temporal','amplitude','cli','regression','output'):
        p.add_argument('--'+name+'-dir',type=Path,required=True)
    p.add_argument('--contact-temporal-dir',type=Path)
    p.add_argument('--vorticity-dir',type=Path)
    p.add_argument('--vorticity-temporal-dir',type=Path)
    a=p.parse_args();out=a.output_dir;out.mkdir(parents=True,exist_ok=True)
    contact=samples(a.contact_dir,False);acoustic=samples(a.acoustic_dir,True)
    if a.vorticity_dir:contact+=samples(a.vorticity_dir,False)
    time=samples(a.temporal_dir,True);amplitude=samples(a.amplitude_dir,True)
    if a.contact_temporal_dir:time+=samples(a.contact_temporal_dir,False)
    if a.vorticity_temporal_dir:time+=samples(a.vorticity_temporal_dir,False)
    tables={'contact':global_errors(contact),'regions':regional(contact),
            'acoustic':global_errors(acoustic),'acoustic_regions':regional(acoustic),
            'reflection':reflections(a.acoustic_dir),'temporal':temporal(time),
            'temporal_reflection':reflections(a.temporal_dir),'amplitude':global_errors(amplitude),
            'amplitude_reflection':reflections(a.amplitude_dir),'cli':read(a.cli_dir/'cli.csv'),
            'negative':read(a.cli_dir/'negative.csv'),'regression':read(a.regression_dir/'equivalence.csv'),
            'vtk':read(a.regression_dir/'vtk.csv')}
    for name,rows in tables.items():write(out/('primitive_euler_phase9_'+name+'.csv'),rows)
    root=Path(__file__).resolve().parents[1];hashes={}
    directories=[a.contact_dir,a.acoustic_dir,a.temporal_dir,a.amplitude_dir,a.cli_dir,
                 a.contact_temporal_dir,a.vorticity_dir,a.vorticity_temporal_dir]
    for directory in filter(None,directories):
        executable=directory/'boundary_audit'
        if executable.exists():hashes[str(executable)]=hashlib.sha256(executable.read_bytes()).hexdigest()
        for path in sorted(directory.rglob('*')):
            if path.is_file() and path.suffix in ('.csv','.cfg','.su2','.log'):
                hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    changed=subprocess.check_output(['git','diff','--name-only'],cwd=root,text=True).splitlines()
    new=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=root,text=True).splitlines()
    for name in changed+new:
        path=root/name
        if path.suffix not in ('.csv','.json','.md') and path.is_file():hashes[name]=hashlib.sha256(path.read_bytes()).hexdigest()
    exe=root/'build/SU2_CFD/src/SU2_CFD';hashes[str(exe)]=hashlib.sha256(exe.read_bytes()).hexdigest()
    manifest=dict(head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        compiler=subprocess.check_output(['clang++','--version'],text=True).strip(),sha256=hashes)
    (out/'primitive_euler_phase9_provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print({name:len(rows) for name,rows in tables.items()})
    print('Temporal passes',sum(r['passes_five_percent'] for r in tables['temporal']),'/',len(tables['temporal']))


if __name__=='__main__':main()
