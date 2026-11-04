#!/usr/bin/env python3
"""Build/run the real Euler/primitive physical-BC contact comparison."""
import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import subprocess
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--build-dir',type=Path,default=Path('build'))
    p.add_argument('--n',type=int,nargs='+',default=[8,16,32,64])
    p.add_argument('--modes',nargs='+',choices=('wall','wall_oblique','symmetry','farfield','total','mass','supersonic','acoustic_far','acoustic_wall','acoustic_total','acoustic_mass','acoustic_outlet','vorticity_far','vorticity_outlet','vorticity_supersonic'),default=['wall','wall_oblique','symmetry','farfield','total','mass','supersonic'])
    p.add_argument('--families',nargs='+',choices=('regular','jitter'),default=['regular','jitter'])
    p.add_argument('--divisors',type=int,nargs='+',default=[1])
    p.add_argument('--amplitude',type=float,default=1e-5)
    p.add_argument('--build-only',action='store_true')
    p.add_argument('--skip-build',action='store_true')
    args=p.parse_args();root=Path(__file__).resolve().parents[1];build=(root/args.build_dir).resolve()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True);exe=out/'boundary_audit'
    if not args.skip_build:
        entry=next(e for e in json.loads((build/'compile_commands.json').read_text()) if e['file'].endswith('CPhysicalGeometry.cpp'))
        command=shlex.split(entry['command']);flags=[];i=1
        while i<len(command):
            if command[i] in ('-o','-MF','-MQ','-c'):i+=2;continue
            if command[i]!='-MD':flags.append(command[i])
            i+=1
        obj=out/'boundary.o'
        subprocess.run([command[0],*flags,'-I'+str(root/'UnitTests/SU2_CFD'),'-c',str(root/'UnitTests/SU2_CFD/primitive_euler_boundary_audit.cpp'),'-o',str(obj)],cwd=build,check=True)
        lines=subprocess.check_output([str(root/'ninja'),'-C',str(build),'-t','commands','SU2_CFD/src/SU2_CFD'],text=True)
        link=shlex.split(lines.splitlines()[-1]);deps=link[link.index('SU2_CFD/src/libSU2core.a'):]
        subprocess.run([command[0],'-o',str(exe),str(obj),*deps],cwd=build,check=True)
    if args.build_only:return
    records=[];executions=[]
    for family in args.families:
        for mode in args.modes:
            for n in args.n:
                for divisor in args.divisors:
                    label=f'{family}_{mode}_{n}_{divisor}';target=out/(label+'.csv');start=time.monotonic()
                    with (out/(label+'.log')).open('w') as log:
                        result=subprocess.run([str(exe),family,str(n),mode,str(divisor),str(target),str(args.amplitude)],env={**os.environ,'OMP_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT)
                    seconds=time.monotonic()-start
                    if result.returncode:raise RuntimeError('Audit failed: '+str(out/(label+'.log')))
                    with target.open() as stream:records.extend(csv.DictReader(stream))
                    executions.append(dict(case=label,seconds=seconds,returncode=result.returncode))
                    print(label,seconds,flush=True)
    for name,rows in (('contact',records),('execution',executions)):
        with (out/(name+'.csv')).open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


if __name__=='__main__':main()
