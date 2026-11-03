#!/usr/bin/env python3
"""Build diagnostic extensions of the unchanged 6B.1/6B.2 harnesses in scratch.

Anchor-checked instrumentation adds exports, not numerical operations. The complete
generated translation unit and input hashes are kept alongside the executable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess


SIGNED_EXPORT = r'''
    // Diagnostic-only export: production residual and already computed defects.
    std::ofstream signedField(std::string(argv[3])+".signed_"+std::to_string(time));
    signedField << std::setprecision(17);
    double sumError=0;
    for (unsigned long i=0;i<grid.xy.size();++i) {
      const double v=SU2_TYPE::GetValue(ng->nodes->GetVolume(i));
      signedField << i << ' ' << grid.xy[i][0] << ' ' << grid.xy[i][1] << ' ' << v << ' ' << grid.layer(i);
      for (const auto* values:{&tauW,&jTauW,&tauU})
        for (auto value:(*values)[i]) signedField << ' ' << value;
      for (unsigned k=0;k<4;++k) {
        const double sum=-(result.faces[i][k]+result.volume[i][k]+result.boundary[i][k])/v+result.source[i][k];
        sumError=std::max(sumError,std::abs(sum-result.rhs[i][k]));
      }
      // Export each primitive contribution with its actual RHS sign and units.
      for (const auto* values:{&result.faces,&result.volume,&result.boundary})
        for (auto value:(*values)[i]) signedField << ' ' << -value/v;
      for (auto value:result.source[i]) signedField << ' ' << value;
      for (auto value:result.rhs[i]) signedField << ' ' << value;
      signedField << '\n';
    }
    Require(sumError<5e-13,"primitive contribution sum differs from productive RHS");
    diagnostic("contribution_sum_error",time,sumError);
'''


def replace_once(source, anchor, replacement):
    if source.count(anchor) != 1:
        raise ValueError('Harness changed; review instrumentation anchor: '+anchor)
    return source.replace(anchor, replacement)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['defect', 'evolution'])
    parser.add_argument('--time', type=float, choices=[0.1, 0.2], default=0.1)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--build-dir', default='build')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = (root / args.build_dir).resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    original = root / ('UnitTests/SU2_CFD/euler_mms_audit.cpp' if args.mode == 'defect'
                       else 'UnitTests/SU2_CFD/euler_mms_evolution.cpp')
    source = original.read_text()
    if args.mode == 'defect':
        source = replace_once(source, '    diagnostic("init_U_error",time,initU);',
                              SIGNED_EXPORT+'\n    diagnostic("init_U_error",time,initU);')
    elif args.time == 0.2:
        source = replace_once(source, 'constexpr double T=0.1;', 'constexpr double T=0.2;')
        source = replace_once(source, 'std::ceil(T*n*n/0.25)', '2*std::ceil(0.1*n*n/0.25)')
    if args.mode == 'evolution':
        source = replace_once(source,
            "    for (const auto* values:{&wp,&up,&wc,&uc}) for (auto value:(*values)[i]) field << ' ' << value;",
            "    for (const auto* values:{&wp,&up,&wc,&uc}) for (auto value:(*values)[i]) field << ' ' << value;\n"
            "    const auto reference=M::Evaluate(grid.xy[i][0],grid.xy[i][1],T).w;\n"
            "    for (auto value:reference) field << ' ' << value;\n"
            "    for (auto value:M::Conservative(reference)) field << ' ' << value;")
    generated = output / 'instrumented.cpp'
    generated.write_text(source)
    entry = next(e for e in json.loads((build / 'compile_commands.json').read_text())
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
    obj, exe = output / 'instrumented.o', output / ('euler_asymptotic_'+args.mode)
    subprocess.run([command[0], *flags, '-I'+str(root), '-I'+str(original.parent),
                    '-c', str(generated), '-o', str(obj)], cwd=build, check=True)
    commands = subprocess.check_output([str(root/'ninja'), '-C', str(build), '-t', 'commands',
                                       'SU2_CFD/src/SU2_CFD'], text=True)
    link = shlex.split(commands.splitlines()[-1])
    dependencies = link[link.index('SU2_CFD/src/libSU2core.a'):]
    subprocess.run([command[0], '-o', str(exe), str(obj), *dependencies], cwd=build, check=True)
    (output/'build_manifest.json').write_text(json.dumps(dict(mode=args.mode, time=args.time,
        original=str(original), original_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
        generated_sha256=hashlib.sha256(generated.read_bytes()).hexdigest(),
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()), indent=2)+'\n')
    print(exe)


if __name__ == '__main__':
    main()
