#!/usr/bin/env python3
"""Build a source-free vortex extension of the unchanged paired FE harness.

Generated source is retained in scratch. All replacements are anchor-checked;
neither productive operators nor previous harnesses are edited.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
from run_euler_asymptotic_study import replace_once

SNAPSHOT = r'''
  const auto snapshot=[&](double t) {
    std::ofstream file(std::string(argv[4])+".snapshot_"+std::to_string(t));
    file << std::setprecision(17);
    for (unsigned long i=0;i<grid.xy.size();++i) {
      file << i << ' ' << SU2_TYPE::GetValue(ng->nodes->GetVolume(i)) << ' ' << grid.xy[i][0] << ' ' << grid.xy[i][1];
      for (const auto* values:{&wp,&up,&wc,&uc}) for (auto value:(*values)[i]) file << ' ' << value;
      const auto reference=M::Evaluate(grid.xy[i][0],grid.xy[i][1],t).w;
      for (auto value:reference) file << ' ' << value;
      for (auto value:M::Conservative(reference)) file << ' ' << value;
      file << '\n';
    }
    if (t==0) {
      std::ofstream cells(std::string(argv[4])+".triangles");
      for (const auto& triangle:grid.triangles) cells << triangle[0] << ' ' << triangle[1] << ' ' << triangle[2] << '\n';
    }
  };
  snapshot(0);
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--build-dir', default='build')
    a = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = (root/a.build_dir).resolve()
    out = a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    old = root/'UnitTests/SU2_CFD/euler_mms_evolution.cpp'
    s = old.read_text()
    s = replace_once(s, '#include "euler_mms_test_support.hpp"', '#include "euler_mms_test_support.hpp"\n#include "euler_isentropic_vortex.hpp"')
    s = replace_once(s, 'namespace M=EulerCommonMMS;', 'namespace M=IsentropicVortex;')
    s = replace_once(s, 'bool IsManufacturedSolution() const override { return true; }', 'bool IsManufacturedSolution() const override { return false; }')
    s = replace_once(s, 'const auto value=M::ConservativeSource(M::Evaluate(SU2_TYPE::GetValue(x[0]),SU2_TYPE::GetValue(x[1]),SU2_TYPE::GetValue(t)));',
        'throw std::runtime_error("No manufactured source is permitted for the vortex");\n    const Vec value{};')
    s = replace_once(s, 'return M::PrimitiveSource(M::Evaluate(q[0],q[1],t));',
        'throw std::runtime_error("No primitive source callback is permitted");\n    return Vec{};')
    s = replace_once(s, 'Grid grid(n,family);',
        'Grid grid(n,family);\n  for (auto& x:grid.xy) for (auto& coordinate:x) coordinate=10*coordinate-5;')
    s = replace_once(s, 'n>=8 && divisor>0 && (family=="regular" || family=="smooth" || family=="jitter")',
        '(n==8 || n==16 || n==32 || n==64) && divisor>0 && (family=="regular" || family=="jitter")')
    s = replace_once(s, 'constexpr double T=0.1;', 'constexpr double T=M::finalTime;')
    s = replace_once(s, 'std::ceil(T*n*n/0.25)', 'std::ceil(T*n*n/(0.02*M::length*M::length))')
    s = replace_once(s, 'Require(exact->sourceCalls==steps*grid.xy.size() && sourceCallsP==steps*grid.xy.size(),"source stage count mismatch");',
        'Require(exact->sourceCalls==0 && sourceCallsP==0,"vortex must have no source calls");')
    s = replace_once(s, 'prim->AdvanceExplicitEuler(pg,pc,stageTime,dt,boundary,source)',
        'prim->AdvanceExplicitEuler(pg,pc,stageTime,dt,boundary)')
    s = replace_once(s, 'q0p=sum(up);q0c=sum(uc);', 'q0p=sum(up);q0c=sum(uc);\n'+SNAPSHOT)
    s = replace_once(s, '    if (step%std::max(1u,steps/4)==0)',
        '    if (step+1==steps/2 || step+1==steps) snapshot(nextTime);\n    if (step%std::max(1u,steps/4)==0)')
    s = replace_once(s, 'if (grid.in(i,g)) {',
        'if (g==5 ? std::hypot(grid.xy[i][0]-M::center0-T,grid.xy[i][1]-M::center0-T)<=2 : grid.in(i,g)) {')
    s = replace_once(s, "<< groupNames[g] << ',' << k", "<< (g==5?\"moving\":groupNames[g]) << ',' << k")
    generated = out/'vortex.cpp'
    generated.write_text(s)
    entry = next(e for e in json.loads((build/'compile_commands.json').read_text())
                 if e['file'].endswith('CPhysicalGeometry.cpp'))
    command = shlex.split(entry['command'])
    flags, i = [], 1
    while i < len(command):
        if command[i] in ('-o', '-MF', '-MQ', '-c'):
            i += 2
            continue
        if command[i] != '-MD':
            flags.append(command[i])
        i += 1
    obj, exe = out/'vortex.o', out/'euler_vortex'
    subprocess.run([command[0], *flags, '-I'+str(root), '-I'+str(old.parent), '-c', str(generated), '-o', str(obj)], cwd=build, check=True)
    commands = subprocess.check_output([str(root/'ninja'), '-C', str(build), '-t', 'commands', 'SU2_CFD/src/SU2_CFD'], text=True)
    link = shlex.split(commands.splitlines()[-1])
    subprocess.run([command[0], '-o', str(exe), str(obj), *link[link.index('SU2_CFD/src/libSU2core.a'):]], cwd=build, check=True)
    (out/'manifest.json').write_text(json.dumps(dict(original_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
        generated_sha256=hashlib.sha256(generated.read_bytes()).hexdigest(), executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()), indent=2)+'\n')
    print(exe)


if __name__ == '__main__':
    main()
