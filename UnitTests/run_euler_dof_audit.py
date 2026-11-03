#!/usr/bin/env python3
"""Build/run the isolated CEulerSolver DOF experiment using native SU2 archives."""
import argparse
import json
from pathlib import Path
import shlex
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', default='build')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--run', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = (root / args.build_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
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
    source = root / 'UnitTests/SU2_CFD/euler_dof_audit.cpp'
    obj, exe = output / 'audit.o', output / 'euler_dof_audit'
    subprocess.run([command[0], *flags, '-I' + str(root), '-c', str(source), '-o', str(obj)],
                   cwd=build, check=True)
    commands = subprocess.check_output([str(root / 'ninja'), '-C', str(build), '-t', 'commands',
                                        'SU2_CFD/src/SU2_CFD'], text=True)
    link = shlex.split(commands.splitlines()[-1])
    dependencies = link[link.index('SU2_CFD/src/libSU2core.a'):]
    subprocess.run([command[0], '-o', str(exe), str(obj), *dependencies], cwd=build, check=True)
    print(exe, flush=True)
    if args.run is not None:
        subprocess.run([str(exe), *args.run], check=True)


if __name__ == '__main__':
    main()
