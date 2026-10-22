#!/usr/bin/env python3
"""Build a focused Catch2 driver using an existing SU2 build; optionally run CLI smoke tests.

--asan instruments the focused tests, factory, primitive classes and variable bases.
It does NOT turn the existing SU2/MPI archives into a fully sanitized build.
"""
import argparse
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', default='build')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--asan', action='store_true')
    parser.add_argument('--cli', action='store_true')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = (root / args.build_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.cli:
        cli_checks(build, output)
        return
    database = json.loads((build / 'compile_commands.json').read_text())
    entry = next(item for item in database if item['file'].endswith('CPhysicalGeometry.cpp'))
    command = shlex.split(entry['command'])
    flags = []
    index = 1
    while index < len(command):
        token = command[index]
        if token in ('-o', '-MF', '-MQ', '-c'):
            index += 2
            continue
        if token != '-MD':
            flags.append(token)
        index += 1
    flags += ['-I' + str(root / 'externals/catch2'), '-I' + str(root)]
    sanitize = ['-fsanitize=address', '-fno-omit-frame-pointer', '-O1', '-g'] if args.asan else []
    if args.asan and sys.platform == 'darwin':
        # Existing native archives were built with libc++ container annotations
        # disabled. Match that setting across template instantiations to avoid
        # mixed-instrumentation poisoning. Heap/stack ASan checks remain enabled;
        # logical vector-capacity checks require a fully instrumented build.
        sanitize.append('-D_LIBCPP_HAS_NO_ASAN')
    runner = output / 'main.cpp'
    runner.write_text('''#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "Common/include/parallelization/mpi_structure.hpp"
int main(int argc, char** argv) {
  SU2_MPI::Init(&argc, &argv);
  int status = Catch::Session().run(argc, argv);
  SU2_MPI::Finalize();
  return status;
}
''')
    sources = [runner] + [root / 'UnitTests/SU2_CFD' / name for name in (
        'primitive_euler_geometry_tests.cpp', 'numerics/primitive_euler_2d_tests.cpp',
        'primitive_euler_solver_tests.cpp')]
    if args.asan:
        sources += [root / 'SU2_CFD/src' / name for name in (
            'solvers/CSolverFactory.cpp', 'solvers/CPrimitiveEulerSolver.cpp',
            'variables/CPrimitiveEulerVariable.cpp', 'variables/CVariable.cpp', 'variables/CScalarVariable.cpp')]
    objects = []
    for number, source in enumerate(sources):
        obj = output / (str(number) + '.o')
        cmd = [command[0]] + flags + sanitize + ['-c', str(source), '-o', str(obj)]
        print(shlex.join(cmd), flush=True)
        subprocess.run(cmd, cwd=build, check=True)
        objects.append(str(obj))
    # Reuse the actual CFD link dependencies; exclude its main() object.
    ninja_commands = subprocess.check_output(
        [str(root / 'ninja'), '-C', str(build), '-t', 'commands', 'SU2_CFD/src/SU2_CFD'], text=True)
    link = shlex.split(ninja_commands.splitlines()[-1])
    dependencies = link[link.index('SU2_CFD/src/libSU2core.a'):]
    cmd = [command[0]] + sanitize + ['-o', str(output / 'test_driver')] + objects + dependencies
    print(shlex.join(cmd), flush=True)
    subprocess.run(cmd, cwd=build, check=True)


def cli_checks(build, output):
    mesh = '''NDIME=2
NELEM=1
5 0 1 2
NPOIN=3
0 0 0
2 0 1
0 1 2
NMARK=1
MARKER_TAG=boundary
MARKER_ELEMS=3
3 0 1
3 1 2
3 2 0
'''
    (output / 'mesh.su2').write_text(mesh)
    common = '''MGLEVEL=0
TIME_DISCRE_FLOW=EULER_EXPLICIT
MESH_FILENAME=mesh.su2
MARKER_FAR=(boundary)
ITER=1
CFL_NUMBER=0.1
CONV_NUM_METHOD_FLOW=ROE
MACH_NUMBER=0.2
MUSCL_FLOW=NO
SLOPE_LIMITER_FLOW=NONE
OUTPUT_FILES=(RESTART)
'''
    checks = [('primitive', 'PRIMITIVE_EULER', '', 1, 'PRIMITIVE_EULER phase 3: storage/factory available'),
              ('euler', 'EULER', '', 0, 'Exit Success (SU2_CFD)'),
              ('harmonic', 'PRIMITIVE_EULER', 'TIME_DOMAIN=YES\nTIME_MARCHING=HARMONIC_BALANCE\nTIME_STEP=0.01\n',
               1, 'Harmonic Balance not yet implemented'),
              ('species', 'PRIMITIVE_EULER', 'KIND_SCALAR_MODEL=SPECIES_TRANSPORT\n',
               1, 'Species transport currently only available')]
    for label, solver, extra, expected, message in checks:
        config = output / (label + '.cfg')
        config.write_text('SOLVER=' + solver + '\n' + common + extra)
        result = subprocess.run([str(build / 'SU2_CFD/src/SU2_CFD'), config.name], cwd=output,
                                capture_output=True, text=True)
        log = result.stdout + result.stderr
        (output / (label + '.log')).write_text(log)
        print(label, 'exit', result.returncode, 'expected', expected, flush=True)
        if result.returncode != expected or message not in log:
            raise RuntimeError('Unexpected CLI result; inspect ' + str(output / (label + '.log')))
        if label == 'primitive' and 'Exit Success (SU2_CFD)' in log:
            raise RuntimeError('Incomplete primitive solver reported success')
        if label == 'euler' and not re.search(r'^\|\s*0\s*\|', log, re.MULTILINE):
            raise RuntimeError('Euler did not report its requested iteration')


if __name__ == '__main__':
    main()
