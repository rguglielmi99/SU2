#!/usr/bin/env python3
"""Run native physical-BC examples with both production solvers and check W."""
import argparse
import csv
import os
from pathlib import Path
import subprocess
import time
import re
import numpy as np
from primitive_euler_boundary_cases import case, NAMES
from run_primitive_euler_driver_e2e import binary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();exe=a.executable.resolve();out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for family in ('regular','jitter'):
        for name in (*NAMES,'temperature_farfield','mass_inlet'):
            original='farfield_subsonic' if name=='temperature_farfield' else ('inlet_outlet_subsonic' if name=='mass_inlet' else name)
            cfg,grid=case(original,family,8)
            if name=='temperature_farfield': cfg=cfg.replace('FREESTREAM_OPTION=DENSITY_FS','FREESTREAM_OPTION=TEMPERATURE_FS')+'FREESTREAM_TEMPERATURE=1\n'
            if name=='mass_inlet':
                cfg=cfg.replace('INLET_TYPE=TOTAL_CONDITIONS','INLET_TYPE=MASS_FLOW')
                cfg=re.sub(r'MARKER_INLET=\(inlet,[^\n]+','MARKER_INLET=(inlet,1,0.4,1,0,0)',cfg)
            solvers=('PRIMITIVE_EULER','EULER','NAVIER_STOKES') if family=='regular' and name=='wall' else ('PRIMITIVE_EULER','EULER')
            for solver in solvers:
                path=out/(family+'_'+name+'_'+solver);path.mkdir(exist_ok=True)
                (path/(original+'.su2')).write_text(grid)
                (path/'case.cfg').write_text(cfg.replace('SOLVER=PRIMITIVE_EULER','SOLVER='+solver))
                start=time.monotonic()
                with (path/'run.log').open('w') as log:
                    result=subprocess.run([str(exe),'case.cfg'],cwd=path,env={**os.environ,'OMP_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT)
                if result.returncode:raise RuntimeError(str(path/'run.log'))
                final=binary(sorted(path.glob('*_state_*.dat'))[-1])
                if solver=='PRIMITIVE_EULER':
                    w=np.column_stack([final[k] for k in ('Density','Velocity_x','Velocity_y','Pressure')])
                else:
                    rho=final['Density'];u=final['Momentum_x']/rho;v=final['Momentum_y']/rho
                    w=np.column_stack((rho,u,v,.4*(final['Energy']-.5*rho*(u*u+v*v))))
                speed=2 if 'supersonic' in name else .4
                error=float(np.max(abs(w-np.array([1,speed,0,1]))))
                if error>2e-13:raise ValueError('Uniform state changed: '+str(path)+' '+str(error))
                rows.append(dict(family=family,case=name,solver=solver,max_state_error=error,min_rho=float(w[:,0].min()),min_p=float(w[:,3].min()),seconds=time.monotonic()-start))
                print(family,name,solver,error,flush=True)
    with (out/'cli.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    negative=[]
    base,grid=case('inlet_outlet_subsonic','regular',8)
    far,fargrid=case('farfield_subsonic','regular',8)
    variants={
        'legacy_far':(far.replace('PRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC','PRIMITIVE_BOUNDARY_MODE=UNIFORM'),fargrid,'native physical freestream'),
        'far_requires_opt_in':(far.replace('PRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC\n',''),fargrid,'LOCAL_CHARACTERISTIC'),
        'us_units':(far.replace('FLUID_MODEL=IDEAL_GAS','FLUID_MODEL=STANDARD_AIR')+'SYSTEM_MEASUREMENTS=US\n',fargrid,'SYSTEM_MEASUREMENTS=SI'),
        'bad_outlet_pressure':(base.replace('MARKER_OUTLET=(outlet,1)','MARKER_OUTLET=(outlet,-1)'),grid,'positive static pressure'),
        'bad_inlet_direction':(base.replace(',1,0,0)',',-1,0,0)'),grid,'direction must point'),
        'subsonic_supersonic_outlet':(base.replace('MARKER_OUTLET=(outlet,1)','MARKER_SUPERSONIC_OUTLET=(outlet)'),grid,'no incoming characteristics'),
        'backflow':(far.replace('MARKER_FAR=(lateral,inlet,outlet)','MARKER_FAR=(lateral,inlet)\nMARKER_OUTLET=(outlet,1)')+'AOA=180\n',fargrid,'backflow'),
        'periodic':(far.replace('MARKER_FAR=(lateral,inlet,outlet)','MARKER_FAR=(lateral)\nMARKER_PERIODIC=(inlet,outlet,0,0,0,0,0,0,10,0,0)'),fargrid,'periodic geometry/residual coupling'),
    }
    for name,(cfg,grid,message) in variants.items():
        path=out/('negative_'+name);path.mkdir(exist_ok=True)
        filename=next(line.split('=',1)[1] for line in cfg.splitlines() if line.startswith('MESH_FILENAME='))
        (path/filename).write_text(grid);(path/'case.cfg').write_text(cfg)
        result=subprocess.run([str(exe),'case.cfg'],cwd=path,env={**os.environ,'OMP_NUM_THREADS':'1'},stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        (path/'run.log').write_text(result.stdout)
        if result.returncode==0 or message not in result.stdout:raise RuntimeError(str(path/'run.log'))
        negative.append(dict(case=name,returncode=result.returncode,expected=message))
    with (out/'negative.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(negative[0]));writer.writeheader();writer.writerows(negative)


if __name__=='__main__':main()
