#!/usr/bin/env python3
"""Bounded phase-15B diagnostic only: 50-step smoke then <=500 local iterations.

Uses the ordinary SU2_CFD driver, phase-15A config, and native output writers.
The archived phase-15A physical run is read only for comparison, never replaced.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from run_primitive_euler_driver_e2e import binary
from run_primitive_euler_naca15 import ROOT, CASES, MESH, configuration, setting, run, snapshots, pvd, table, write, vtu


def config(steps, cfl):
    text=configuration('primitive', steps)
    text=re.sub(r'(?m)^MAX_TIME=.*\n', '', text)
    for key,value in [('PRIMITIVE_LOCAL_PSEUDO_TIME','YES'),('UNST_CFL_NUMBER',0),('CFL_NUMBER',cfl),
        ('HISTORY_OUTPUT','(ITER,PRIMITIVE_RHS,PRIMITIVE_STATE,AERO_COEFF,PRIMITIVE_STEADY_AUDIT,LOCAL_PSEUDO_TIME)')]:
        text=setting(text,key,value)
    return text


def check(folder):
    history=table(folder/'history.csv')
    records=snapshots(folder,'primitive')
    pvd(folder,records)
    assert 'Cur_Time' not in history[0] and 'Time_Step' not in history[0]
    for row in history:
        assert row['Min_rho']>0 and row['Min_p']>0 and row['MAX_MACH']<1 and row['MAX_TRACE_MACH']<1
        if row['Time_Iter']>0:
            assert row['DT_PSEUDO_MIN']>0 and row['DT_PSEUDO_MAX']>=row['DT_PSEUDO_MIN']
            np.testing.assert_allclose(row['DT_PSEUDO_RATIO'],row['DT_PSEUDO_MAX']/row['DT_PSEUDO_MIN'],rtol=1e-14)
            assert row['DT_PSEUDO_MIN']<=row['DT_PSEUDO_TE']<=row['DT_PSEUDO_MAX']
    for record in records:
        _,fields,_=vtu(folder/record['sequence_file'])
        assert 'Physical_Time' not in fields
        np.testing.assert_array_equal(fields['Pseudo_Iteration'],record['iteration'])
        if record['iteration']>0: assert np.min(fields['Local_Pseudo_Step'])>0
    last=history[-1]
    _,fields,_=vtu(folder/records[-1]['sequence_file'])
    steps=fields['Local_Pseudo_Step']
    write(folder/'step_distribution.csv',[dict(percentile=p,step=float(np.percentile(steps,p))) for p in [0,10,25,50,75,90,100]])
    execution=json.loads((folder/'execution.json').read_text())
    return dict(run=folder.name,iterations=int(last['Time_Iter']),seconds=execution['seconds'],
        seconds_per_iteration=execution['seconds']/last['Time_Iter'],
        residual_initial=history[0]['RHS_CURRENT'],residual_final=last['RHS_CURRENT'],
        residual_ratio=last['RHS_CURRENT']/history[0]['RHS_CURRENT'],
        interior_initial=history[0]['RHS_INTERIOR'],interior_final=last['RHS_INTERIOR'],
        te_initial=history[0]['RHS_TE'],te_final=last['RHS_TE'],
        dt_min=last['DT_PSEUDO_MIN'],dt_max=last['DT_PSEUDO_MAX'],dt_ratio=last['DT_PSEUDO_RATIO'],dt_te=last['DT_PSEUDO_TE'],
        effective_cfl=last['EFFECTIVE_CFL'],rejected_iterations=last['PSEUDO_REJECTED'],
        rejected_attempts=last['PSEUDO_REJECTED_ATTEMPTS'],max_mach=max(r['MAX_MACH'] for r in history),
        max_trace_mach=max(r['MAX_TRACE_MACH'] for r in history),min_rho=min(r['Min_rho'] for r in history),
        min_p=min(r['Min_p'] for r in history),CL=last['CL'],CD=last['CD'])


def physical_regression(exe,out,seconds):
    baseline=CASES/'results/smoke_primitive'
    folder=out/'physical_smoke5'
    run(exe,folder,(baseline/'case.cfg').read_text(),seconds)
    snapshots(folder,'primitive')
    differences=[]
    for original in sorted(baseline.glob('primitive_state*.dat')):
        old,new=binary(original),binary(folder/original.name)
        assert old.keys()==new.keys()
        differences.extend(float(np.max(np.abs(old[key]-new[key]))) for key in old)
    assert differences and max(differences)<2e-13
    (out/'physical_regression.json').write_text(json.dumps(dict(max_state_difference=max(differences),
        steps=5,baseline=str(baseline)),indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable',type=Path,default=ROOT/'build/SU2_CFD/src/SU2_CFD')
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--iterations',type=int,default=300)
    parser.add_argument('--cfl',type=float,default=.05)
    parser.add_argument('--seconds-per-run',type=float,default=120)
    parser.add_argument('--postprocess-only',action='store_true')
    args=parser.parse_args()
    if not 200<=args.iterations<=500 or not np.isfinite(args.cfl) or args.cfl<=0 or not np.isfinite(args.seconds_per_run) or args.seconds_per_run<=0:
        parser.error('Use 200..500 iterations, a finite positive CFL and a positive time budget')
    out=args.output_dir.resolve()
    checksum=hashlib.sha256(MESH.read_bytes()).hexdigest()
    assert checksum==(CASES/'results/mesh.sha256').read_text().strip()
    if not args.postprocess_only:
        out.mkdir(parents=True,exist_ok=False)
        physical_regression(args.executable.resolve(),out,args.seconds_per_run)
        run(args.executable.resolve(),out/'smoke50',config(50,args.cfl),args.seconds_per_run)
        smoke=check(out/'smoke50')
        cfl=args.cfl
        if smoke['residual_ratio']>1:
            cfl*=.5
            run(args.executable.resolve(),out/'smoke50_lower_cfl',config(50,cfl),args.seconds_per_run)
            reduced=check(out/'smoke50_lower_cfl')
            if reduced['residual_ratio']>1:
                print(json.dumps(dict(smoke=smoke,reduced=reduced,status='Residual grows: no extended run'),indent=2))
                write(out/'summary.csv',[smoke,reduced])
                return
        run(args.executable.resolve(),out/'diagnostic',config(args.iterations,cfl),args.seconds_per_run)
    rows=[check(folder) for folder in sorted(out.iterdir()) if folder.is_dir() and folder.name!='physical_smoke5' and (folder/'history.csv').exists()]
    write(out/'summary.csv',rows)
    baseline=table(CASES/'results/primitive/history.csv')
    matching=[]
    local=table(out/'diagnostic/history.csv') if (out/'diagnostic/history.csv').exists() else []
    for count in [0,50,200,300,500,3000]:
        for name,history in [('global_FE',baseline),('local_pseudo',local)]:
            matches=[r for r in history if int(r['Time_Iter'])==count]
            if matches:
                r=matches[0];matching.append(dict(method=name,iteration=count,rhs=r['RHS_CURRENT'],rhs_te=r['RHS_TE'],
                    CL=r['CL'],CD=r['CD'],min_rho=r['Min_rho'],min_p=r['Min_p']))
    write(out/'comparison_global_fe.csv',matching)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    main()
