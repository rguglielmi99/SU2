#!/usr/bin/env python3
"""Phase 15D: staged, bounded warm-start/TE experiments. No native solver runs.

Pilots first; final runs require an explicit initialization, TE policy and CFL.
All runs start from freestream or the original native field, never a silent
continuation. Budget counts all solver subprocesses across invocations.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from convert_primitive_euler_initial import convert
from run_primitive_euler_driver_e2e import binary
from run_primitive_euler_local_pseudo15b import config
from run_primitive_euler_naca15 import ROOT,CASES,MESH,setting,write,vtu
from run_primitive_euler_steady15c import history,checkpoints,REGIONS,compare_cp,marker_ids


def run(exe,out,label,initial,policy,cfl,steps,budget):
    executions=[json.loads(p.read_text()) for p in out.glob('run_*/execution.json')]
    remaining=budget-sum(r['seconds'] for r in executions)
    if remaining<15: raise ValueError('Computational budget exhausted')
    folder=out/('run_'+label); folder.mkdir(exist_ok=False)
    text=config(steps,cfl)
    for key,value in [('PRIMITIVE_TE_POLICY',policy),('OUTPUT_FILES','(PARAVIEW_MULTIBLOCK,PARAVIEW,RESTART)'),
                      ('OUTPUT_WRT_FREQ','(1000,300,300)'),('SCREEN_WRT_FREQ_TIME',1000)]:
        text=setting(text,key,value)
    if initial=='native':
        text=setting(text,'PRIMITIVE_INITIALIZATION','FILE')
        text=setting(text,'PRIMITIVE_INITIAL_FILE','../native_initial.dat')
    (folder/'case.cfg').write_text(text)
    start=time.monotonic(); reason='iteration_cap'
    with (folder/'run.log').open('w') as log:
        p=subprocess.Popen([str(exe),'case.cfg'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,
                           env={**os.environ,'OMP_NUM_THREADS':'1'})
        while p.poll() is None:
            rows=history(folder/'history.csv')
            if time.monotonic()-start>remaining-8: reason='time_budget'
            if rows:
                r=rows[-1]
                if not all(np.isfinite(v) for v in r.values()) or r['Min_rho']<=0 or r['Min_p']<=0 or r['MAX_TRACE_MACH']>=1 or r['MAX_MACH']>=1:
                    reason='physical_guard'
                if len(rows)>800:
                    means=np.mean(np.array([r['RHS_CURRENT'] for r in rows[-800:]]).reshape(4,200),axis=1)
                    if np.all(means[1:]>1.05*means[:-1]): reason='persistent_growth'
                (folder/'progress.json').write_text(json.dumps(dict(iteration=r['Time_Iter'],rhs=r['RHS_CURRENT'],
                    te=r['RHS_TE'],CL=r['CL'],CD=r['CD'],seconds=time.monotonic()-start)))
            if reason!='iteration_cap':
                p.terminate()
                try:p.wait(timeout=5)
                except subprocess.TimeoutExpired:p.kill();p.wait();reason+=':forced_kill'
                break
            time.sleep(.5)
    record=dict(label=label,initial=initial,policy=policy,cfl=cfl,requested_steps=steps,
                seconds=time.monotonic()-start,returncode=p.returncode,stop_reason=reason)
    (folder/'execution.json').write_text(json.dumps(record,indent=2))
    print(json.dumps(record),flush=True)
    if p.returncode!=0: raise RuntimeError('Solver failed; inspect '+str(folder/'run.log'))
    checkpoints(folder)


def projection(out,folder):
    raw=binary(out/'native_initial.dat'); loaded=binary(folder/'primitive_state_00000.dat')
    ids=np.arange(len(loaded['Density']))
    assert len(ids)==5233
    for key in ['x','y']: np.testing.assert_array_equal(loaded[key],raw[key])
    initial=np.column_stack([raw[k] for k in ['Density','Velocity_x','Velocity_y','Pressure']])[ids]
    actual=np.column_stack([loaded[k] for k in ['Density','Velocity_x','Velocity_y','Pressure']])
    np.testing.assert_array_equal(actual[:,[0,3]],initial[:,[0,3]])
    wall=np.isin(ids,list(marker_ids()['airfoil'])); tip=ids==199
    rows=[]
    for name,mask in [('ordinary_wall',wall&~tip),('te',tip),('other',~wall)]:
        diff=actual[mask]-initial[mask]
        rows.append(dict(region=name,nodes=int(np.sum(mask)),max_state_correction=float(np.max(abs(diff))),
            max_velocity_correction=float(np.max(np.linalg.norm(diff[:,1:3],axis=1))),
            max_rho_correction=float(np.max(abs(diff[:,0]))),max_p_correction=float(np.max(abs(diff[:,3])))))
        if name=='other': np.testing.assert_array_equal(diff,0)
    write(folder/'initial_projection.csv',rows)
    # Full precision difference from the imported native field, at initialization
    # and the last available restart. This is a drift diagnostic, not an error
    # against an exact solution or an equivalence claim between formulations.
    last=max(folder.glob('primitive_state_*.dat'),key=lambda p:int(p.stem.rsplit('_',1)[1]))
    point_fields=vtu(folder/'primitive_00000.vtu')[1]
    volumes=np.empty(len(ids));volumes[point_fields['Global_ID'].astype(int)]=point_fields['Dual_Volume']
    differences=[]
    for path in [folder/'primitive_state_00000.dat',last]:
        field=binary(path)
        for key in ['x','y']:np.testing.assert_array_equal(field[key],raw[key])
        for key in ['Density','Velocity_x','Velocity_y','Pressure']:
            delta=field[key]-raw[key]
            differences.append(dict(iteration=int(path.stem.rsplit('_',1)[1]),component=key,
                max_abs_difference=float(np.max(abs(delta))),
                volume_weighted_rms_difference=float(np.sqrt(np.sum(volumes*delta**2)/np.sum(volumes)))))
    write(folder/'field_difference_native.csv',differences)


def localize(folder,step,out,label):
    # Data-only analysis of saved phase-15C/15D fields. Not an operator change.
    targets=list(folder.glob(f'primitive_{step:05d}.vtu'))+list(folder.glob(f'primitive_{step:06d}.vtu'))
    if not targets:
        import xml.etree.ElementTree as ET
        for path in folder.glob('primitive_*.vtm'):
            if int(path.stem.rsplit('_',1)[1])==step:
                targets=[folder/node.attrib['file'] for node in ET.parse(path).getroot().findall('.//DataSet') if node.attrib['name']=='Internal']
    xyz,f,_=vtu(targets[0]); ids=f['Global_ID'].astype(int); wall=np.isin(ids,list(marker_ids()['airfoil'])); far=np.isin(ids,list(marker_ids()['farfield']))
    norm=np.sum(np.column_stack([f[k] for k in ['RhsRho','RhsU','RhsV','RhsP']]).astype(float)**2,axis=1)
    weighted=f['Dual_Volume'].astype(float)*norm
    interior=~(wall|far); nearte=np.linalg.norm(xyz[:,:2]-[1,0],axis=1)<.1
    nearle=np.linalg.norm(xyz[:,:2],axis=1)<.1
    rows=[]
    for name,mask in [('interior_TE_r01',interior&nearte),('interior_LE_r01',interior&nearle),
                      ('interior_other',interior&~(nearte|nearle)),('wall',wall),('farfield',far)]:
        rows.append(dict(case=label,iteration=step,region=name,nodes=int(sum(mask)),
            global_square_fraction=float(sum(weighted[mask])/sum(weighted)),
            rhs=float(np.sqrt(sum(weighted[mask])/sum(f['Dual_Volume'][mask])))))
    write(out/(label+'_localization.csv'),rows)
    top=np.argsort(-weighted)[:30]
    write(out/(label+'_largest_residual_nodes.csv'),[dict(GlobalIndex=int(ids[i]),x=float(xyz[i,0]),y=float(xyz[i,1]),
        rhs=float(np.sqrt(norm[i])),global_square_fraction=float(weighted[i]/sum(weighted))) for i in top])


def process(out):
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/primitive-phase15-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    baseline=CASES/'results15c/cfl_0p05'; histories={'freestream_WEAK_reference':history(baseline/'history.csv')[:301]}
    summaries=[]; check=[]
    for folder in sorted(out.glob('run_*')):
        if not (folder/'execution.json').exists():continue
        rows,saved=checkpoints(folder); first,last=rows[0],rows[-1]
        rec=json.loads((folder/'execution.json').read_text()); histories[folder.name]=rows
        window=rows[-min(200,len(rows)):]
        s=dict(**rec,iterations=int(last['Time_Iter']),saved_iteration=saved,
            rhs_initial=first['RHS_CURRENT'],rhs_final=last['RHS_CURRENT'],
            decades=float(np.log10(first['RHS_CURRENT']/last['RHS_CURRENT'])),
            freestream_decades=float(np.log10(.04702700066616893/last['RHS_CURRENT'])),
            interior_initial=first['RHS_INTERIOR'],interior_final=last['RHS_INTERIOR'],
            te_initial=first['RHS_TE'],te_final=last['RHS_TE'],wall_final=last['RHS_WALL_STRONG'],farfield_final=last['RHS_FARFIELD'],
            CL=last['CL'],CD=last['CD'],CL_range=float(np.ptp([r['CL'] for r in window])),CD_range=float(np.ptp([r['CD'] for r in window])),
            min_rho=min(r['Min_rho'] for r in rows),min_p=min(r['Min_p'] for r in rows),
            max_mach=max(r['MAX_MACH'] for r in rows),max_trace_mach=max(r['MAX_TRACE_MACH'] for r in rows),
            rejected_attempts=last['PSEUDO_REJECTED_ATTEMPTS'],dt_min=last['DT_PSEUDO_MIN'],dt_max=last['DT_PSEUDO_MAX'])
        s['decades_per_second']=s['decades']/s['seconds']
        # Warm starts must also have a small ABSOLUTE residual, not merely a small initial update.
        s['stationarity_gate']=bool(s['rhs_final']<.04702700066616893*1e-3 and
            s['interior_final']<.026036125094*1e-3 and s['te_final']<1e-3 and
            s['CL_range']<1e-5 and s['CD_range']<1e-5 and saved==s['iterations'])
        summaries.append(s)
        selected={0,300,1000,2000,3000,5000,6000,9000,int(last['Time_Iter'])}
        for r in rows:
            if int(r['Time_Iter']) not in selected:continue
            check.append(dict(run=folder.name,**r,**{g+'_normalized':r[k]/first[k] if first[k]>1e-12 else '' for g,k in REGIONS.items()}))
        if rec['initial']=='native':projection(out,folder)
        localize(folder,0,out,folder.name+'_initial')
        localize(folder,saved,out,folder.name+'_final')
        if s['stationarity_gate']: compare_cp(out,folder,saved)
    write(out/'summary.csv',summaries);write(out/'checkpoints.csv',check)
    write(out/'baseline_reference.csv',[dict(source=str(baseline),**r) for r in [histories['freestream_WEAK_reference'][0],histories['freestream_WEAK_reference'][-1]]])
    for field in ['RHS_CURRENT','RHS_INTERIOR','RHS_WALL_STRONG','RHS_TE','CL','CD']:
        fig,ax=plt.subplots()
        for label,rows in histories.items():
            (ax.semilogy if field.startswith('RHS') else ax.plot)([r['Time_Iter'] for r in rows],[r[field] for r in rows],label=label)
        ax.set(xlabel='Pseudo-iteration',ylabel=field); ax.grid();ax.legend(fontsize=6);fig.tight_layout();fig.savefig(out/(field+'.png'));plt.close(fig)
    fig,ax=plt.subplots()
    for label,rows in histories.items():
        ax.semilogy([r['Time_Iter'] for r in rows],[r['RHS_CURRENT']/rows[0]['RHS_CURRENT'] for r in rows],label=label)
    ax.set(xlabel='Pseudo-iteration',ylabel='Current spatial RHS / own initial RHS')
    ax.grid();ax.legend(fontsize=6);fig.tight_layout();fig.savefig(out/'RHS_normalized.png');plt.close(fig)
    for name,step in [('cfl_0p05',3000),('cfl_0p1',3000),('cfl_0p2',6000)]:
        localize(CASES/'results15c'/name,step,out,'phase15c_'+name)
    (out/'cp_comparison_status.txt').write_text('Stationarity gate: '+str(any(s['stationarity_gate'] for s in summaries))+'\n')
    (out/'budget.json').write_text(json.dumps(dict(solver_seconds=sum(s['seconds'] for s in summaries),
        note='All 15D pilots/finals included; reused baseline and offline postprocessing excluded'),indent=2))
    print(json.dumps(summaries,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=CASES/'results15d')
    p.add_argument('--executable',type=Path,default=ROOT/'build/SU2_CFD/src/SU2_CFD')
    p.add_argument('--stage',choices=['pilots','final','postprocess'],required=True)
    p.add_argument('--label',default='final');p.add_argument('--initial',choices=['native','freestream'],default='native')
    p.add_argument('--policy',choices=['WEAK','PROJECTED'],default='PROJECTED')
    p.add_argument('--cfl',type=float,default=.05);p.add_argument('--steps',type=int,default=9000)
    p.add_argument('--budget-seconds',type=float,default=480)
    a=p.parse_args();out=a.output_dir.resolve()
    if not 0<a.cfl<=.2 or not 0<a.budget_seconds<=480 or a.steps<1:p.error('Require CFL in (0,.2], budget <=480 s, positive steps')
    if a.stage=='pilots':
        out.mkdir(parents=True,exist_ok=False)
        native_config=CASES/'results/native/case.cfg'
        options={line.split('=',1)[0].strip():line.split('=',1)[1].strip().split('%')[0].strip()
                 for line in native_config.read_text().splitlines() if '=' in line and not line.lstrip().startswith('%')}
        assert options['SOLVER']=='EULER' and options['REF_DIMENSIONALIZATION']=='DIMENSIONAL'
        for key,value in [('GAMMA_VALUE',1.4),('GAS_CONSTANT',1),('MACH_NUMBER',.2),('AOA',0),
                          ('FREESTREAM_DENSITY',1),('FREESTREAM_PRESSURE',1),('REF_AREA',1)]:
            assert float(options[key])==value,(key,options[key])
        report=convert(CASES/'results/native/native_state_000599.dat',MESH,out/'native_initial.dat')
        report['native_config_sha256']=hashlib.sha256(native_config.read_bytes()).hexdigest()
        (out/'conversion.json').write_text(json.dumps(report,indent=2))
        for initial,policy in [('freestream','PROJECTED'),('native','WEAK'),('native','PROJECTED')]:
            run(a.executable.resolve(),out,'pilot_'+initial+'_'+policy,initial,policy,.05,300,a.budget_seconds)
    elif a.stage=='final':
        run(a.executable.resolve(),out,a.label,a.initial,a.policy,a.cfl,a.steps,a.budget_seconds)
    process(out)


if __name__=='__main__':main()
