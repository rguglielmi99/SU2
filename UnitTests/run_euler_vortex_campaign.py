#!/usr/bin/env python3
"""Progressive native vortex campaign; refuses smooth and N>64 by construction."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from euler_vortex_output import write_vtu, write_pvd


def validate_case(path, n):
    cells = np.loadtxt(str(path)+'.triangles', dtype=int)
    if cells.shape != (2*n*n, 3):
        raise ValueError('Wrong triangle count')
    for t in (0., .25, .5):
        a = np.loadtxt(str(path)+f'.snapshot_{t:.6f}')
        if a.shape != ((n+1)**2, 28) or not np.all(np.isfinite(a)):
            raise ValueError('Invalid snapshot')
        np.testing.assert_array_equal(a[:, 0], np.arange((n+1)**2))
        if abs(a[:, 1].sum()-100)>1e-11:
            raise ValueError('Wrong dual volume sum')
        if t == 0 and max(np.max(abs(a[:, 4:8]-a[:, 20:24])), np.max(abs(a[:, 16:20]-a[:, 24:28]))) > 2e-14:
            raise ValueError('Wrong analytic initialization')
    with Path(str(path)+'.checks').open() as f:
        check = {r['kind']: float(r['value']) for r in csv.DictReader(f)}
    for key in ('source_calls_primitive', 'source_calls_conservative', 'native_fallbacks', 'contractions'):
        if check[key] != 0:
            raise ValueError('Source/fallback/contraction in vortex: '+key)
    if check['T'] != .5 or check['Dirichlet_error'] > 2e-14 or check['theta_min'] != 1:
        raise ValueError('Time/BC/positivity check failed')
    return check


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    action = p.add_mutually_exclusive_group(required=True)
    action.add_argument('--pilot', action='store_true')
    action.add_argument('--run', action='store_true')
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--divisors', nargs='+', type=int, default=[1, 2, 4])
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    if any(d not in (1, 2, 4, 8, 16) for d in a.divisors):
        p.error('Unsupported divisor')
    out = a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    sha = hashlib.sha256(a.executable.read_bytes()).hexdigest()
    manifest = out/'pilot_passed.json'
    if manifest.exists() and json.loads(manifest.read_text())['executable_sha256'] != sha:
        raise ValueError('Executable changed: choose a fresh output directory; do not mix trajectories')

    def run(task):
        family, n, d = task
        path = out/f'{family}_{n}_{d}.csv'
        if a.resume and path.exists():
            validate_case(path, n)
            return f'{task}: validated existing'
        start = time.monotonic()
        with path.with_suffix('.log').open('w') as log:
            subprocess.run([str(a.executable.resolve()), family, str(n), str(d), str(path)],
                check=True, stdout=log, stderr=subprocess.STDOUT, env={**os.environ, 'OMP_NUM_THREADS': '1'})
        validate_case(path, n)
        return f'{task}: passed in {time.monotonic()-start:.2f}s'

    if a.pilot:
        print(run(('regular', 16, 1)), flush=True)
        path = out/'regular_16_1.csv'
        cells = np.loadtxt(str(path)+'.triangles', dtype=int)
        entries = []
        for i, t in enumerate((0., .25, .5)):
            field = np.loadtxt(str(path)+f'.snapshot_{t:.6f}')
            name = f'pilot_{i}.vtu'
            write_vtu(out/name, field, cells, t)
            entries.append((t, name))
        write_pvd(out/'pilot.pvd', entries)
        manifest.write_text(json.dumps(dict(executable_sha256=sha, pilot='regular N=16', writer='XML roundtrip passed'), indent=2)+'\n')
        print('Pilot and VTU/PVD validation passed; campaign may proceed.')
        return
    if not manifest.exists() or json.loads(manifest.read_text())['executable_sha256'] != sha:
        raise ValueError('Run the pilot with this executable before the campaign')
    tasks = [(f, n, d) for f in ('regular', 'jitter') for n in (8, 16, 32, 64) for d in a.divisors]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for message in pool.map(run, tasks):
            print(message, flush=True)


if __name__ == '__main__':
    main()
