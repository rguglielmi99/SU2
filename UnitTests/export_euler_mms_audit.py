#!/usr/bin/env python3
"""Run and collect the phase-6B.1 semidiscrete experiment (no time stepping)."""
import argparse
import concurrent.futures
import csv
import math
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--collect-only', action='store_true')
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    tasks = [(f, n) for f in ('regular', 'smooth', 'jitter') for n in (8, 16, 32, 64)]

    def run(task):
        family, n = task
        dest = output / f'{family}_{n}.csv'
        with dest.with_suffix('.log').open('w') as log:
            subprocess.run([str(args.executable.resolve()), family, str(n), str(dest)],
                           stdout=log, stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, 'OMP_NUM_THREADS': '1'})
        return f'{family} N={n}: passed'

    if not args.collect_only:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for status in pool.map(run, tasks):
                print(status, flush=True)
    rows = []
    for f, n in tasks:
        with (output / f'{f}_{n}.csv').open() as stream:
            rows.extend(csv.DictReader(stream))
    diagnostics = [r for r in rows if r['group'] == 'check']
    errors = [r for r in rows if r['group'] != 'check']
    previous = {}
    for r in errors:
        key = tuple(r[k] for k in ('kind', 'family', 'time', 'group', 'component'))
        for norm in ('rms', 'l1', 'l2', 'linf'):
            value = float(r[norm])
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'Invalid norm: {r}')
            old = previous.get((key, norm))
            r['p_'+norm] = (math.log(old[1]/value)/math.log(int(r['N'])/old[0])
                            if r['kind'] != 'Dirichlet_constraint' and old and old[1] > 0 and value > 0 else '')
            previous[key, norm] = int(r['N']), value
    for filename, data in [('errors.csv', errors), ('checks.csv', diagnostics)]:
        with (output / filename).open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f'Exported {len(errors)} norm rows and {len(diagnostics)} checks', flush=True)


if __name__ == '__main__':
    main()
