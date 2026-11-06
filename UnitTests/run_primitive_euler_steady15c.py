#!/usr/bin/env python3
"""Budgeted phase-15C campaign, unchanged spatial solver, independent freestream starts.

Only the primitive solver is executed. Native phase-15A results are read only
if the explicit stationarity gate is met. No restart continuation is assumed.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET

import numpy as np
from run_primitive_euler_local_pseudo15b import config
from run_primitive_euler_naca15 import ROOT, CASES, MESH, setting, table, write, vtu, pvd
from run_primitive_euler_driver_e2e import binary

REGIONS = {'global': 'RHS_CURRENT', 'interior': 'RHS_INTERIOR', 'wall': 'RHS_WALL_STRONG',
           'te': 'RHS_TE', 'farfield': 'RHS_FARFIELD'}


def history(path):
    if not path.exists(): return []
    text = path.read_text()
    lines = text.splitlines()
    if not text.endswith('\n'): lines = lines[:-1]
    if len(lines)<2: return []
    reader = csv.DictReader(lines, skipinitialspace=True)
    reader.fieldnames = [key.strip() for key in reader.fieldnames]
    return [{key:float(value) for key,value in row.items()} for row in reader]


def run(exe, folder, cfl, steps, seconds, total_deadline):
    folder.mkdir(exist_ok=False)
    text = config(steps,cfl)
    # Native writers: VTM at thousands, standalone VTU at multiples of 300.
    # They use different filenames; selected PVD checkpoints are deduplicated.
    text = setting(text,'OUTPUT_FILES','(PARAVIEW_MULTIBLOCK,PARAVIEW,RESTART)')
    text = setting(text,'OUTPUT_WRT_FREQ','(1000,300,1000)')
    text = setting(text,'SCREEN_WRT_FREQ_TIME',1000)
    (folder/'case.cfg').write_text(text)
    start = time.monotonic()
    deadline = min(start+seconds,total_deadline-12)
    cause = 'iteration_cap'
    with (folder/'run.log').open('w') as log:
        proc = subprocess.Popen([str(exe),'case.cfg'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,
                                env={**os.environ,'OMP_NUM_THREADS':'1'})
        while proc.poll() is None:
            rows = history(folder/'history.csv')
            if time.monotonic()>=deadline: cause='time_budget'
            if rows:
                last = rows[-1]
                if last['Min_rho']<=0 or last['Min_p']<=0 or last['MAX_MACH']>=1 or last['MAX_TRACE_MACH']>=1:
                    cause='physical_guard'
                if len(rows)>=801:
                    values=np.array([r['RHS_CURRENT'] for r in rows[-800:]])
                    blocks=np.mean(values.reshape(4,200),axis=1)
                    if np.all(blocks[1:]>1.05*blocks[:-1]): cause='persistent_residual_growth'
                (folder/'progress.json').write_text(json.dumps(dict(iteration=last['Time_Iter'],
                    rhs=last['RHS_CURRENT'],rhs_te=last['RHS_TE'],CD=last['CD'],seconds=time.monotonic()-start)))
            if cause!='iteration_cap':
                proc.terminate()  # SU2 STOP handler writes the last accepted state.
                try: proc.wait(timeout=10)
                except subprocess.TimeoutExpired: proc.kill(); proc.wait(); cause+=':forced_kill'
                break
            time.sleep(.5)
    result=dict(cfl=cfl,requested_iterations=steps,seconds=time.monotonic()-start,
                returncode=proc.returncode,stop_reason=cause,initial_state='freestream')
    (folder/'execution.json').write_text(json.dumps(result,indent=2))
    return result


def marker_ids():
    lines=MESH.read_text().splitlines()
    result={}
    for i,line in enumerate(lines):
        if line.startswith('MARKER_TAG='):
            tag=line.split('=')[1].strip(); count=int(lines[i+1].split('=')[1])
            result[tag]={int(k) for edge in lines[i+2:i+2+count] for k in edge.split()[1:3]}
    return result


def checkpoints(folder):
    rows=table(folder/'history.csv'); by_step={int(row['Time_Iter']):row for row in rows}
    targets={}
    for path in sorted(folder.glob('primitive_*.vtu')):
        targets[int(path.stem.rsplit('_',1)[1])]=path
    for path in sorted(folder.glob('primitive_*.vtm')):
        nodes=ET.parse(path).getroot().findall('.//DataSet')
        for node in nodes:
            assert (folder/node.attrib['file']).is_file()
            if node.attrib['name']=='Internal': targets[int(path.stem.rsplit('_',1)[1])]=folder/node.attrib['file']
    assert targets
    masks=marker_ids(); records=[]; regional=[]
    selected={0,300,1000,2000,3000,5000,6000,max(targets)}
    for step,path in sorted(targets.items()):
        xyz,fields,cells=vtu(path)
        assert len(xyz)==5233 and cells==10216
        ids=fields['Global_ID'].astype(int)
        np.testing.assert_array_equal(np.sort(ids),np.arange(5233))
        assert 'Physical_Time' not in fields
        np.testing.assert_array_equal(fields['Pseudo_Iteration'],step)
        if step>0: assert np.min(fields['Local_Pseudo_Step'])>0
        assert np.min(fields['Density'])>0 and np.min(fields['Pressure'])>0
        mach=np.linalg.norm(fields['Velocity'][:,:2],axis=1)/np.sqrt(1.4*fields['Pressure']/fields['Density'])
        assert np.max(mach)<1
        np.testing.assert_allclose(mach,fields['Mach'],rtol=5e-7,atol=5e-8)
        volume=fields['Dual_Volume'].astype(float)
        norm=np.sum(np.column_stack([fields[k] for k in ['RhsRho','RhsU','RhsV','RhsP']]).astype(float)**2,axis=1)
        wall=np.isin(ids,list(masks['airfoil'])); far=np.isin(ids,list(masks['farfield'])); tip=ids==199
        groups={'global':np.ones(len(ids),bool),'interior':~(wall|far),'wall':wall&~tip,'te':tip,'farfield':far}
        np.testing.assert_array_equal(sum(groups[g].astype(int) for g in ['interior','wall','te','farfield']),1)
        total=float(np.sum(volume*norm)); partition=0
        for group,mask in groups.items():
            weighted=float(np.sum(volume[mask]*norm[mask])); weight=float(np.sum(volume[mask]))
            rms=np.sqrt(weighted/weight)
            np.testing.assert_allclose(rms,by_step[step][REGIONS[group]],rtol=3e-6,atol=1e-13)
            if group!='global': partition+=weighted
            initial=rows[0][REGIONS[group]]
            regional.append(dict(iteration=step,region=group,nodes=int(np.sum(mask)),volume=weight,rhs=rms,
                normalized=rms/initial if initial>1e-12 else '',
                global_square_fraction=weighted/total if total else 0))
        np.testing.assert_allclose(partition,total,rtol=1e-14)
        if step in selected: records.append(dict(iteration=step,file=str(path.relative_to(folder)),sequence_file=str(path.relative_to(folder))))
    pvd(folder,records)
    write(folder/'regional_checkpoints.csv',regional)
    last_fields=vtu(targets[max(targets)])[1]
    write(folder/'step_distribution.csv',[dict(percentile=p,step=float(np.percentile(last_fields['Local_Pseudo_Step'],p))) for p in [0,10,25,50,75,90,100]])
    return rows,max(targets)


def process(out):
    os.environ.setdefault('MPLCONFIGDIR','/private/tmp/primitive-phase15-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summaries=[]; checkpoints_table=[]; histories={}
    for folder in sorted(out.glob('cfl_*')):
        if not (folder/'history.csv').exists(): continue
        rows,saved=checkpoints(folder); first,last=rows[0],rows[-1]; histories[folder.name]=rows
        execution=json.loads((folder/'execution.json').read_text()); window=rows[-min(200,len(rows)):]
        summary=dict(run=folder.name,**execution,iterations=int(last['Time_Iter']),saved_iteration=saved,
            rhs_initial=first['RHS_CURRENT'],rhs_final=last['RHS_CURRENT'],decades=np.log10(first['RHS_CURRENT']/last['RHS_CURRENT']),
            interior=last['RHS_INTERIOR'],wall=last['RHS_WALL_STRONG'],te=last['RHS_TE'],farfield=last['RHS_FARFIELD'],
            interior_decades=np.log10(first['RHS_INTERIOR']/last['RHS_INTERIOR']),
            CL=last['CL'],CD=last['CD'],CD_window_range=max(r['CD'] for r in window)-min(r['CD'] for r in window),
            CL_window_range=max(r['CL'] for r in window)-min(r['CL'] for r in window),
            min_rho=min(r['Min_rho'] for r in rows),min_p=min(r['Min_p'] for r in rows),
            max_mach=max(r['MAX_MACH'] for r in rows),max_trace_mach=max(r['MAX_TRACE_MACH'] for r in rows),
            dt_min=last['DT_PSEUDO_MIN'],dt_max=last['DT_PSEUDO_MAX'],dt_te=last['DT_PSEUDO_TE'],
            effective_cfl=last['EFFECTIVE_CFL'],rejected_attempts=last['PSEUDO_REJECTED_ATTEMPTS'],
            decades_per_second=np.log10(first['RHS_CURRENT']/last['RHS_CURRENT'])/execution['seconds'])
        # A conservative reporting gate, not a solver convergence modification.
        summary['stationarity_gate']=bool(summary['returncode']==0 and saved==summary['iterations'] and
            summary['decades']>=3 and summary['interior_decades']>=3 and summary['te']<1e-3 and
            summary['te']/first['RHS_TE']<=1e-3 and summary['CD_window_range']<1e-5 and summary['CL_window_range']<1e-5)
        summaries.append(summary)
        requested={0,300,1000,2000,3000,5000,6000,int(last['Time_Iter'])}
        for r in rows:
            if int(r['Time_Iter']) not in requested: continue
            row=dict(run=folder.name,iteration=int(r['Time_Iter']))
            for k in ['RHS_CURRENT','RHS_INTERIOR','RHS_WALL_STRONG','RHS_TE','RHS_FARFIELD','CL','CD','Min_rho','Min_p',
                      'MAX_MACH','MAX_TRACE_MACH','DT_PSEUDO_MIN','DT_PSEUDO_MAX','DT_PSEUDO_TE','EFFECTIVE_CFL','PSEUDO_REJECTED_ATTEMPTS']:
                row[k]=r[k]
            for region,key in REGIONS.items():
                row[region+'_normalized']=r[key]/first[key] if first[key]>1e-12 else ''
            te_file=folder/f'surface_primitive_te_{row["iteration"]}.csv'
            te=table(te_file)[0] if te_file.exists() else {}
            for k in ['rho','u','v','p']: row['TE_'+k]=te.get(k,'')
            checkpoints_table.append(row)
        # Windows quantify plateaus/oscillations without assuming exponential decay.
        windows=[]
        for begin in range(1,len(rows),200):
            block=rows[begin:begin+200]
            windows.append(dict(first=block[0]['Time_Iter'],last=block[-1]['Time_Iter'],
                rhs_mean=np.mean([r['RHS_CURRENT'] for r in block]),te_mean=np.mean([r['RHS_TE'] for r in block]),
                wall_mean=np.mean([r['RHS_WALL_STRONG'] for r in block]),interior_mean=np.mean([r['RHS_INTERIOR'] for r in block]),
                CD_range=np.ptp([r['CD'] for r in block]),CL_range=np.ptp([r['CL'] for r in block])))
        write(folder/'windows.csv',windows)
    write(out/'summary.csv',summaries); write(out/'checkpoints.csv',checkpoints_table)
    for field,name in [('RHS_CURRENT','residuals'),('RHS_TE','trailing_edge'),('CL','lift'),('CD','drag')]:
        fig,ax=plt.subplots()
        for label,rows in histories.items():
            x=[r['Time_Iter'] for r in rows]; y=np.array([r[field] for r in rows])
            if field=='RHS_CURRENT': y/=y[0]
            (ax.semilogy if field.startswith('RHS') else ax.plot)(x,y,label=label)
        ax.set(xlabel='Pseudo-iteration',ylabel=field+(' / initial' if field=='RHS_CURRENT' else ''))
        ax.grid(); ax.legend(); fig.tight_layout(); fig.savefig(out/(name+'.png')); plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax in axes.flat: ax.set_visible(False)
    for ax,(label,rows) in zip(axes.flat,histories.items()):
        ax.set_visible(True)
        for group,key in REGIONS.items():
            if rows[0][key]>1e-12:
                ax.semilogy([r['Time_Iter'] for r in rows],[r[key]/rows[0][key] for r in rows],label=group)
        ax.set(title=label,xlabel='Pseudo-iteration',ylabel='Regional RHS / initial'); ax.grid(); ax.legend()
    fig.tight_layout(); fig.savefig(out/'regions.png'); plt.close(fig)
    # No quantitative Cp comparison is generated unless stationarity is established.
    (out/'cp_comparison_status.txt').write_text('Stationarity gate: '+str(any(s['stationarity_gate'] for s in summaries))+'\n')
    for summary in summaries:
        if summary['stationarity_gate']: compare_cp(out,out/summary['run'],summary['saved_iteration'])
    print(json.dumps(summaries,indent=2))


def compare_cp(out,folder,step):
    import matplotlib.pyplot as plt
    with (folder/f'surface_primitive_pressure_{step}.csv').open() as stream:
        primitive=np.array([[float(r[k]) for k in ['x','y','Cp']] for r in csv.DictReader(stream)])
    native_folder=CASES/'results/native'
    path=sorted(native_folder.glob('native_*.vtm'))[-1]
    node=next(n for n in ET.parse(path).getroot().findall('.//DataSet')
              if n.attrib['name']=='airfoil' and Path(n.attrib['file']).parts[0]==path.stem)
    xyz,fields,_=vtu(native_folder/node.attrib['file'])
    native=np.column_stack([xyz[:,:2],(fields['Pressure'].astype(float)-1)/.028])
    rows=[]; fig,ax=plt.subplots()
    for side,sign in [('upper',1),('lower',-1)]:
        a=primitive[sign*primitive[:,1]>0]; b=native[sign*native[:,1]>=0]
        a=a[np.argsort(a[:,0])]; b=b[np.argsort(b[:,0])]
        assert np.all(np.diff(b[:,0])>0)
        a=a[(a[:,0]>=b[0,0])&(a[:,0]<=b[-1,0])]
        ref=np.interp(a[:,0],b[:,0],b[:,2])
        for value,exact in zip(a,ref): rows.append(dict(side=side,x=value[0],primitive_cp=value[2],native_cp=exact,difference=value[2]-exact))
        ax.plot(a[:,0],a[:,2],label='primitive '+side); ax.plot(b[:,0],b[:,2],'--',label='native '+side)
    write(folder/'cp_comparison.csv',rows)
    difference=np.array([r['difference'] for r in rows])
    write(folder/'cp_metrics.csv',[dict(points=len(rows),rms=np.sqrt(np.mean(difference**2)),maximum=np.max(abs(difference)))])
    ax.invert_yaxis(); ax.set(xlabel='x/c',ylabel='Cp'); ax.legend(); ax.grid();fig.tight_layout();fig.savefig(folder/'cp.png');plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,default=ROOT/'build/SU2_CFD/src/SU2_CFD')
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--budget-seconds',type=float,default=480)
    p.add_argument('--postprocess-only',action='store_true')
    a=p.parse_args(); out=a.output_dir.resolve()
    if not np.isfinite(a.budget_seconds) or a.budget_seconds<60 or a.budget_seconds>480:
        p.error('Use an explicit total budget between 60 and 480 seconds')
    assert hashlib.sha256(MESH.read_bytes()).hexdigest()==(CASES/'results/mesh.sha256').read_text().strip()
    if not a.postprocess_only:
        out.mkdir(parents=True,exist_ok=False); start=time.monotonic(); deadline=start+a.budget_seconds
        preflight=out/'preflight50'
        result=run(a.executable.resolve(),preflight,.05,50,15,deadline)
        assert result['returncode']==0
        checkpoints(preflight)
        original=CASES/'results15b/smoke50'
        diffs=[]
        for path in original.glob('primitive_state*.dat'):
            target=preflight/path.name
            if target.exists():
                old,new=binary(path),binary(target)
                diffs.extend(np.max(abs(old[k]-new[k])) for k in old)
        assert diffs and max(diffs)==0
        (out/'output_regression.json').write_text(json.dumps(dict(max_state_difference=float(max(diffs))),indent=2))
        for cfl,steps,allocation in [(.05,3000,140),(.1,3000,130),(.2,6000,480)]:
            remaining=deadline-time.monotonic()
            if remaining<25: break
            folder=out/('cfl_'+str(cfl).replace('.','p'))
            result=run(a.executable.resolve(),folder,cfl,steps,min(allocation,remaining-12),deadline)
            print(folder.name,result,flush=True)
            if result['stop_reason']=='physical_guard': break
        (out/'campaign_execution.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,budget=a.budget_seconds),indent=2))
    process(out)


if __name__=='__main__': main()
