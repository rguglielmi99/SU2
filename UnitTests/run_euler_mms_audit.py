#!/usr/bin/env python3
"""Build isolated common-MMS analytic tests and a semidiscrete SU2 harness."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', default='build')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--target', choices=['both', 'analytic', 'harness'], default='both')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = (root / args.build_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    database = json.loads((build / 'compile_commands.json').read_text())
    command = shlex.split(next(e['command'] for e in database if e['file'].endswith('CPhysicalGeometry.cpp')))
    flags, i = [], 1
    while i < len(command):
        if command[i] in ('-o', '-MF', '-MQ', '-c'):
            i += 2
            continue
        if command[i] != '-MD':
            flags.append(command[i])
        i += 1
    flags += ['-I' + str(root), '-I' + str(root / 'externals/catch2')]
    link = shlex.split(subprocess.check_output(
        [str(root / 'ninja'), '-C', str(build), '-t', 'commands', 'SU2_CFD/src/SU2_CFD'], text=True).splitlines()[-1])
    dependencies = link[link.index('SU2_CFD/src/libSU2core.a'):]
    runner = output / 'analytic_main.cpp'
    runner.write_text('#define CATCH_CONFIG_MAIN\n#include "catch.hpp"\n')
    for label, sources in [('euler_mms_analytic', [runner, root / 'UnitTests/SU2_CFD/euler_mms_analytic_tests.cpp']),
                           ('euler_mms_audit', [root / 'UnitTests/SU2_CFD/euler_mms_audit.cpp'])]:
        if args.target == 'analytic' and label != 'euler_mms_analytic':
            continue
        if args.target == 'harness' and label != 'euler_mms_audit':
            continue
        objects = []
        for number, source in enumerate(sources):
            obj = output / f'{label}_{number}.o'
            subprocess.run([command[0], *flags, '-c', str(source), '-o', str(obj)], cwd=build, check=True)
            objects.append(str(obj))
        subprocess.run([command[0], '-o', str(output / label), *objects, *dependencies], cwd=build, check=True)
        print(output / label, flush=True)


if __name__ == '__main__':
    main()
