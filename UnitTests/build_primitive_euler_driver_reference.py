#!/usr/bin/env python3
"""Add geometry/first-step exports to the existing phase-7 reference harness.

This is a test reference only. The driver under test is the ordinary SU2_CFD.
"""
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys
from run_euler_asymptotic_study import replace_once

GEOMETRY=r'''
  {
    std::ofstream file(std::string(argv[4])+".geometry.csv");
    file<<std::setprecision(17)<<"kind,first,second,x,y,nx,ny,volume\n";
    const auto& metrics=prim->GetGeometryMetrics();
    for (unsigned long i=0;i<ng->GetnPoint();++i) {
      const auto& cell=metrics.GetCell(i);
      file<<"node,"<<ng->nodes->GetGlobalIndex(i)<<",0,"<<ng->nodes->GetCoord(i,0)<<','<<ng->nodes->GetCoord(i,1)
          <<','<<cell.offset[0]<<','<<cell.offset[1]<<','<<cell.volume<<'\n';
    }
    for (unsigned long e=0;e<metrics.GetnEdge();++e) {
      const auto& f=metrics.GetFace(e);
      file<<"edge,"<<ng->nodes->GetGlobalIndex(f.first)<<','<<ng->nodes->GetGlobalIndex(f.second)<<','
          <<f.point[0]<<','<<f.point[1]<<','<<f.normal[0]<<','<<f.normal[1]<<",0\n";
    }
    for (const auto& f:metrics.GetBoundaryFaces())
      file<<"boundary,"<<ng->nodes->GetGlobalIndex(f.node)<<','<<f.marker<<','<<f.point[0]<<','<<f.point[1]
          <<','<<f.normal[0]<<','<<f.normal[1]<<",0\n";
  }
'''


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--build-dir',default='build')
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1];out=a.output_dir.resolve();build=(root/a.build_dir).resolve()
    subprocess.run([sys.executable,str(root/'UnitTests/run_euler_vortex.py'),'--output-dir',str(out),'--build-dir',a.build_dir],check=True)
    src=out/'vortex.cpp'
    text=replace_once(src.read_text(),'  q0p=sum(up);q0c=sum(uc);','  q0p=sum(up);q0c=sum(uc);\n'+GEOMETRY)
    text=replace_once(text,'if (step+1==steps/2 || step+1==steps) snapshot(nextTime);',
                      'if (step+1==1 || step+1==steps/2 || step+1==steps) snapshot(nextTime);')
    instrumented=out/'reference.cpp';instrumented.write_text(text)
    entry=next(e for e in json.loads((build/'compile_commands.json').read_text()) if e['file'].endswith('CPhysicalGeometry.cpp'))
    command=shlex.split(entry['command']);flags=[];i=1
    while i<len(command):
        if command[i] in ('-o','-MF','-MQ','-c'):i+=2;continue
        if command[i]!='-MD':flags.append(command[i])
        i+=1
    obj=out/'reference.o';exe=out/'reference_vortex'
    subprocess.run([command[0],*flags,'-I'+str(root),'-I'+str(root/'UnitTests/SU2_CFD'),'-c',str(instrumented),'-o',str(obj)],cwd=build,check=True)
    lines=subprocess.check_output([str(root/'ninja'),'-C',str(build),'-t','commands','SU2_CFD/src/SU2_CFD'],text=True)
    link=shlex.split(lines.splitlines()[-1]);dependencies=link[link.index('SU2_CFD/src/libSU2core.a'):]
    subprocess.run([command[0],'-o',str(exe),str(obj),*dependencies],cwd=build,check=True)
    print(exe)


if __name__=='__main__':main()
