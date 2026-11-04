#!/usr/bin/env python3
"""Package the native-driver campaign without rerunning or changing its data."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--modes-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = {}

    def record(path):
        records[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()

    for directory, names in (
        (args.input_dir, ('execution', 'uniform', 'equivalence', 'errors',
                          'negative', 'geometry_reference', 'vtk')),
        (args.modes_dir, ('modes', 'rhs_history', 'first_step', 'uniform_rhs')),
    ):
        for name in names:
            source = directory / (name + '.csv')
            record(source)
            shutil.copyfile(source, args.output_dir / ('primitive_euler_phase8_' + name + '.csv'))

    for pattern in ('*.cfg', '*.su2', '*.dat', '*geometry*.csv', '*history*.csv'):
        for directory in (args.input_dir, args.modes_dir):
            for path in sorted(directory.rglob(pattern)):
                record(path)
    for pattern in ('Common/include/verification/*.hpp',
                    'SU2_CFD/*/iteration/CPrimitiveEulerIteration.*',
                    'SU2_CFD/*/output/CPrimitiveEulerOutput.*',
                    'UnitTests/*primitive_euler_driver*.py',
                    'TestCases/primitive_euler_examples/**/*.cfg',
                    'TestCases/primitive_euler_examples/**/*.su2'):
        for path in sorted(root.glob(pattern)):
            record(path)
    # Include every tracked source change, not just newly introduced classes.
    changed = subprocess.check_output(['git', 'diff', '--name-only'], cwd=root, text=True)
    for name in changed.splitlines():
        path = root / name
        if path.is_file():
            record(path)
    record(root / 'build/SU2_CFD/src/SU2_CFD')
    manifest = {
        'starting_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        'branch': subprocess.check_output(['git', 'branch', '--show-current'], cwd=root, text=True).strip(),
        'compiler': subprocess.check_output(['clang++', '--version'], text=True).strip(),
        'precision': {'native_restart': 'Float64 coordinates and W', 'native_vtu': 'Float32'},
        'sha256': dict(sorted(records.items())),
    }
    target = args.output_dir / 'primitive_euler_phase8_provenance.json'
    target.write_text(json.dumps(manifest, indent=2) + '\n')
    print('Packaged 11 CSV files and', len(records), 'SHA256 entries in', args.output_dir)


if __name__ == '__main__':
    main()
