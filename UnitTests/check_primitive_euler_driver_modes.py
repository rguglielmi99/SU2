#!/usr/bin/env python3
"""Additional native CLI checks: CFL, history semantics, MPI and other solvers."""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
from primitive_euler_driver_cases import mesh,config
from run_primitive_euler_driver_e2e import run,binary,write,reference


def history(path):
    with (path/'history.csv').open() as f:
        reader=csv.DictReader(f,skipinitialspace=True);reader.fieldnames=[k.strip() for k in reader.fieldnames]
        return list(reader)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',required=True,type=Path)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--reference-dir',type=Path)
    a=p.parse_args();exe=a.executable.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    base=config('mesh.su2',8,True);results=[];uniform_rhs=[]
    variants={
        'fixed_no_max': '\n'.join(l for l in base.splitlines() if not l.startswith('MAX_TIME='))+'\n',
        'fixed_no_count': '\n'.join(l for l in base.splitlines() if not l.startswith('TIME_ITER='))+'\n',
        'fixed_shorter_max':base.replace('MAX_TIME=0.080000000000000002','MAX_TIME=0.04'),
        'CFL':base.replace('TIME_STEP=0.01','TIME_STEP=0').replace('UNST_CFL_NUMBER=0','UNST_CFL_NUMBER=0.5').replace('TIME_ITER=8','TIME_ITER=100'),
        'CFL_step_cap':'\n'.join(l for l in base.replace('TIME_STEP=0.01','TIME_STEP=0').replace(
            'UNST_CFL_NUMBER=0','UNST_CFL_NUMBER=0.5').replace('TIME_ITER=8','TIME_ITER=3').splitlines() if not l.startswith('MAX_TIME='))+'\n',
        'far_exterior':base.replace('MARKER_CUSTOM=(boundary)','MARKER_FAR=(boundary)'),
        'green_gauss':base.replace('NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES','NUM_METHOD_GRAD_RECON=GREEN_GAUSS'),
        'least_squares':base.replace('NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES','NUM_METHOD_GRAD_RECON=LEAST_SQUARES'),
        'barth':base.replace('SLOPE_LIMITER_FLOW=NONE','SLOPE_LIMITER_FLOW=BARTH_JESPERSEN'),
        'venkat':base.replace('SLOPE_LIMITER_FLOW=NONE','SLOPE_LIMITER_FLOW=VENKATAKRISHNAN')}
    for name,text in variants.items():
        path=out/name;path.mkdir(exist_ok=True);(path/'mesh.su2').write_text(mesh('regular',8))
        code,seconds=run(exe,path,text)
        if code:raise ValueError(f'{name} failed; inspect {path}/run.log')
        h=history(path);tf=float(h[-1]['Cur_Time']);expected=(3*float(h[1]['Time_Step']) if name=='CFL_step_cap' else .04 if name=='fixed_shorter_max' else .08)
        np.testing.assert_allclose(tf,expected,rtol=0,atol=64*np.finfo(float).eps)
        if name.startswith('CFL'):
            if len(h)>=101:raise ValueError('CFL failed to reach MAX_TIME within step cap')
            if name=='CFL_step_cap' and len(h)!=4:raise ValueError('CFL step cap not respected')
            np.testing.assert_allclose(sum(float(r['Time_Step']) for r in h),tf,rtol=0,atol=1e-14)
        else:
            assert len(h)==(5 if name=='fixed_shorter_max' else 9)
            assert all(float(r['Time_Step'])==.01 for r in h[1:])
        maximum=max(abs(float(r['L2_rhs_'+k])) for r in h for k in ('RHO','U','V','P'))
        bound=64*np.finfo(float).eps*3/(10/8)
        if maximum>bound:raise ValueError('Uniform RHS exceeds a roundoff-scale geometric bound')
        uniform_rhs.append(dict(case=name,max_uniform_rhs=maximum,roundoff_bound=bound))
        results.append(dict(case=name,returncode=code,steps=len(h)-1,final_time=tf,seconds=seconds))
    write(out/'uniform_rhs.csv',uniform_rhs)
    path=out/'positivity_rejection';path.mkdir(exist_ok=True);(path/'mesh.su2').write_text(mesh('regular',8))
    code,seconds=run(exe,path,config('mesh.su2',8,False,steps=1,dt=100))
    if code==0 or 'Forward Euler rejected' not in (path/'run.log').read_text():
        raise ValueError('A non-admissible primitive FE step was not rejected')
    if (path/'primitive_state_00001.dat').exists():raise ValueError('Rejected step produced an accepted checkpoint')
    results.append(dict(case='positivity_rejection',returncode=code,steps=0,final_time=0,seconds=seconds))
    # Check the post-step history against FE state differences, with a roundoff-derived bound.
    path=out/'rhs_history';path.mkdir(exist_ok=True);(path/'mesh.su2').write_text(mesh('regular',16))
    text=config('mesh.su2',16,False,steps=8,dt=.5/512).replace('OUTPUT_WRT_FREQ=4','OUTPUT_WRT_FREQ=1')
    code,seconds=run(exe,path,text)
    if code:raise ValueError('RHS history case failed')
    h=history(path)
    with (path/'geometry.csv').open() as f:volume={int(r['first']):float(r['volume']) for r in csv.DictReader(f) if r['kind']=='node'}
    v=np.array([volume[i] for i in range(289)])
    states=[]
    for i in range(9):
        b=binary(path/f'primitive_state_{i:05d}.dat');states.append(np.column_stack([b[k] for k in ('Density','Velocity_x','Velocity_y','Pressure')]))
    rhs=[]
    for i in range(1,9):
        dt=float(h[i]['Time_Step']);rate=(states[i]-states[i-1])/dt
        norm=np.sqrt(np.sum(v[:,None]*rate*rate,axis=0)/v.sum())
        bound=16*np.finfo(float).eps*np.max(abs(states[i-1]))/dt
        for k,key in enumerate(('RHO','U','V','P')):
            difference=abs(norm[k]-float(h[i]['L2_rhs_'+key]))
            if difference>bound:raise ValueError('History does not describe the accepted FE RHS')
            rhs.append(dict(step=i,component=k,history=float(h[i]['L2_rhs_'+key]),from_state_difference=norm[k],absolute_difference=difference,roundoff_bound=bound))
    write(out/'rhs_history.csv',rhs)
    if a.reference_dir:
        file=next((a.reference_dir/'reference_regular_16').glob('reference.csv.snapshot_0.000977'))
        ref=np.loadtxt(file)
        error=np.max(abs(ref[:,4:8]-states[1]))
        if error>2e-13:raise ValueError('First driver step differs from reference')
        results.append(dict(case='first_step_reference',returncode=0,steps=1,final_time=.5/512,seconds=0))
        write(out/'first_step.csv',[dict(max_state_difference=error,tolerance=2e-13)])
    # Existing conservative solver paths retain their own factories, iteration and output.
    common='\n'.join(l for l in base.splitlines() if not l.startswith('PRIMITIVE_'))+'\n'
    common=common.replace('MARKER_CUSTOM=(boundary)','MARKER_FAR=(boundary)').replace(
        'HISTORY_OUTPUT=(ITER,TIME_DOMAIN,PRIMITIVE_RHS,PRIMITIVE_STATE)','HISTORY_OUTPUT=(ITER,TIME_DOMAIN,RMS_RES)')
    for kind in ('EULER','NAVIER_STOKES'):
        path=out/kind;path.mkdir(exist_ok=True);(path/'mesh.su2').write_text(mesh('regular',8))
        text=common.replace('SOLVER=PRIMITIVE_EULER','SOLVER='+kind)
        if kind=='NAVIER_STOKES':text+='\nVISCOSITY_MODEL=CONSTANT_VISCOSITY\nMU_CONSTANT=1e-5\n'
        code,seconds=run(exe,path,text)
        if code:raise ValueError(f'Compatibility case {kind} failed')
        results.append(dict(case=kind,returncode=code,steps=8,final_time=.08,seconds=seconds))
    mpiexec=shutil.which('mpiexec')
    if mpiexec:
        path=out/'MPI_negative';path.mkdir(exist_ok=True);(path/'mesh.su2').write_text(mesh('regular',8));(path/'case.cfg').write_text(base)
        with (path/'run.log').open('w') as f:
            proc=subprocess.run([mpiexec,'-n','2',str(exe),'case.cfg'],cwd=path,stdout=f,stderr=subprocess.STDOUT,
                env={**os.environ,'OMP_NUM_THREADS':'1'},timeout=60)
        if proc.returncode==0 or 'serial mode' not in (path/'run.log').read_text():raise ValueError('MPI was not explicitly rejected')
        results.append(dict(case='MPI_2_rejected',returncode=proc.returncode,steps=0,final_time=0,seconds=0))
    write(out/'modes.csv',results)
    print('Time modes, RHS history, first step, EULER/NS compatibility and MPI checks passed.')


if __name__=='__main__':main()
