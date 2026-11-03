#!/usr/bin/env python3
"""Run phase-6A.1 and export errors, moments, temporal and initialization checks.

All output goes to --output-dir (use a scratch directory). No solver is modified.
"""
import argparse
import concurrent.futures
import csv
import math
import os
from pathlib import Path
import subprocess
import time

GROUPS = ('all', 'boundary', 'ring1', 'ring2', 'deep3', 'fixed', 'interior')


def belongs(point, n, group):
    i, j = point[0] % (n+1), point[0] // (n+1)
    layer = min(i, j, n-i, n-j)
    return (True, layer == 0, layer == 1, layer == 2, layer >= 3,
            0.25 <= point[2] <= 0.75 and 0.25 <= point[3] <= 0.75, layer >= 1)[GROUPS.index(group)]


def field(path):
    return [[int(parts[0]), *map(float, parts[1:])]
            for line in path.read_text().splitlines() if (parts := line.split())]


def difference(a, b, n, group, k):
    samples = [(x[1], x[4+k]-y[4+k]) for x, y in zip(a, b) if belongs(x, n, group)]
    return (math.sqrt(sum(e*e for _, e in samples)/len(samples)),
            math.sqrt(sum(v*e*e for v, e in samples)/sum(v for v, _ in samples)),
            max(abs(e) for _, e in samples))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--sizes', type=int, nargs='+', default=[8, 16, 32, 64])
    parser.add_argument('--families', nargs='+', default=['regular', 'smooth', 'jitter'])
    parser.add_argument('--collect-only', action='store_true')
    parser.add_argument('--extra-records', type=Path)
    parser.add_argument('--benchmark', choices=['sinusoidal', 'affine'], default='sinusoidal')
    parser.add_argument('--geometry-only', action='store_true')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    tasks = [(f, n, init, d) for f in args.families for n in args.sizes
             for init in ('point', 'mean') for d in (1, 2, 4)]
    if args.geometry_only:
        tasks = [(f, n, 'point', 1) for f in args.families for n in args.sizes]

    def name(f, n, init, d):
        return output / f'{f}_{n}_{init}_{d}.csv'

    def run(task):
        path = name(*task)
        start = time.monotonic()
        with path.with_suffix('.log').open('w') as log:
            command = [str(args.executable.resolve()), *map(str, task), str(path)]
            if args.benchmark == 'affine':
                command.append('affine')
            elif args.geometry_only:
                command.append('geometry')
            subprocess.run(command,
                           stdout=log, stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, 'OMP_NUM_THREADS': '1'})
        return f'{task}: {time.monotonic()-start:.2f}s'

    if not args.collect_only:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for message in pool.map(run, tasks):
                print(message, flush=True)

    records = []
    for task in tasks:
        with name(*task).open() as stream:
            records.extend(csv.DictReader(stream))
    if args.extra_records:
        with args.extra_records.open() as stream:
            records.extend(csv.DictReader(stream))
    records.sort(key=lambda r: (r['family'], int(r['N']), r['initial'], int(r['divisor'])))
    keys = list(records[0])
    with (output / 'all_records.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(records)
    selected = [r for r in records if r['divisor'] == ('1' if args.geometry_only else '4')]
    previous = {}
    for r in selected:
        key = tuple(r[k] for k in ('kind', 'family', 'initial', 'group', 'component'))
        for norm in ('rms', 'l2', 'linf'):
            old = previous.get((key, norm))
            new = float(r[norm])
            r['p_'+norm] = (math.log(old[1]/new)/math.log(int(r['N'])/old[0])
                            if old and old[1] > 0 and new > 0 else '')
            previous[key, norm] = (int(r['N']), new)
    with (output / 'summary.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys+['p_rms', 'p_l2', 'p_linf'])
        writer.writeheader()
        writer.writerows(selected)

    if args.geometry_only:
        print(f'Exported {len(tasks)} geometry checks to {output}', flush=True)
        return

    temporal = []
    for family in args.families:
        for n in args.sizes:
            for init in ('point', 'mean'):
                fields = {d: field(Path(str(name(family, n, init, d))+'.field')) for d in (1, 2, 4)}
                for group in GROUPS:
                    for k in range(4):
                        first = difference(fields[1], fields[2], n, group, k)
                        second = difference(fields[2], fields[4], n, group, k)
                        reference = {r['kind']: r for r in selected if r['family'] == family
                                     and int(r['N']) == n and r['initial'] == init
                                     and r['group'] == group and int(r['component']) == k}
                        for index, norm in enumerate(('rms', 'l2', 'linf')):
                            e = float(reference['Epoint'][norm])
                            d = float(reference['D'][norm])
                            temporal.append(dict(family=family, N=n, initial=init, group=group,
                                component=k, norm=norm, delta12=first[index], delta24=second[index],
                                p_time=math.log2(first[index]/second[index]) if first[index]*second[index] > 0 else '',
                                fraction_Epoint=second[index]/e if e > 0 else '',
                                fraction_D=second[index]/d if d > 0 else ''))
    with (output / 'temporal.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(temporal[0]))
        writer.writeheader()
        writer.writerows(temporal)
    with (output / 'initialization.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['family', 'N', 'group', 'component', 'rms', 'l2', 'linf'])
        for family in args.families:
            for n in args.sizes:
                a = field(Path(str(name(family, n, 'point', 4))+'.field'))
                b = field(Path(str(name(family, n, 'mean', 4))+'.field'))
                for group in GROUPS:
                    for k in range(4):
                        writer.writerow([family, n, group, k, *difference(a, b, n, group, k)])
    print(f'Exported {len(tasks)} runs to {output}', flush=True)


if __name__ == '__main__':
    main()
