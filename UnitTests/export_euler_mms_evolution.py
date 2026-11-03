#!/usr/bin/env python3
"""Run paired Forward Euler MMS trajectories and test temporal independence."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from export_euler_dof_audit import GROUPS, belongs, field

PRIMARY = ('W_primitive', 'U_primitive', 'W_conservative', 'U_conservative')
KINDS = (*PRIMARY, 'difference_W', 'difference_U')
NORMS = ('rms', 'l1', 'l2', 'linf')
NOISE = 1e-12  # Unit scales are one; used only to classify unresolved roundoff-level errors.


def values(row, kind):
    if kind in PRIMARY:
        offset = 4+4*PRIMARY.index(kind)
        return row[offset:offset+4]
    return [row[4+k]-row[12+k] for k in range(4)] if kind == 'difference_W' else [row[8+k]-row[16+k] for k in range(4)]


def delta_norm(a, b, n, group, kind, component):
    samples = [(x[1], values(x, kind)[component]-values(y, kind)[component])
               for x, y in zip(a, b) if belongs(x, n, group)]
    weight = sum(v for v, _ in samples)
    return (math.sqrt(sum(e*e for _, e in samples)/len(samples)),
            sum(v*abs(e) for v, e in samples)/weight,
            math.sqrt(sum(v*e*e for v, e in samples)/weight), max(abs(e) for _, e in samples))


def write(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--jobs', default=2, type=int)
    parser.add_argument('--sizes', nargs='+', type=int, default=[8, 16, 32, 64])
    parser.add_argument('--families', nargs='+', default=['regular', 'smooth', 'jitter'])
    parser.add_argument('--divisors', nargs='+', type=int, default=[1, 2, 4])
    parser.add_argument('--collect-only', action='store_true')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    executable_hash = hashlib.sha256(args.executable.resolve().read_bytes()).hexdigest()
    manifest = output / 'manifest.json'
    if manifest.exists():
        saved = json.loads(manifest.read_text())
        if saved['executable_sha256'] != executable_hash:
            raise ValueError('Executable changed: use a new output directory; do not mix trajectories.')
    else:
        manifest.write_text(json.dumps(dict(executable=str(args.executable.resolve()),
            executable_sha256=executable_hash, created_during='collection' if args.collect_only else 'execution',
            temporal_fraction_limit=0.05, order_change_limit=0.1, roundoff_floor=NOISE), indent=2)+'\n')
    sizes, divisors = sorted(set(args.sizes)), sorted(set(args.divisors))
    if len(divisors) < 3 or any(b != 2*a for a, b in zip(divisors, divisors[1:])):
        parser.error('At least three consecutively doubled divisors are required.')

    def name(f, n, d):
        return output / f'{f}_{n}_{d}.csv'

    tasks = [(f, n, d) for f in args.families for n in sizes for d in divisors]

    def run(task):
        path = name(*task)
        if args.resume and path.exists() and Path(str(path)+'.field').exists() and Path(str(path)+'.checks').exists():
            return f'{task}: existing completed output'
        with path.with_suffix('.log').open('w') as log:
            subprocess.run([str(args.executable.resolve()), *map(str, task), str(path)],
                           stdout=log, stderr=subprocess.STDOUT, check=True,
                           env={**os.environ, 'OMP_NUM_THREADS': '1'})
        return f'{task}: passed'

    if not args.collect_only:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for result in pool.map(run, tasks):
                print(result, flush=True)
    records, checks = [], []
    for task in tasks:
        with name(*task).open() as stream:
            records.extend(csv.DictReader(stream))
        with Path(str(name(*task))+'.checks').open() as stream:
            checks.extend(csv.DictReader(stream))
    previous = {}
    for r in records:
        key = tuple(r[k] for k in ('kind', 'family', 'divisor', 'group', 'component'))
        constrained = r['kind'].endswith('conservative') and r['group'] == 'boundary'
        for norm in NORMS:
            e = float(r[norm])
            if not math.isfinite(e) or e < 0:
                raise ValueError(r)
            old = previous.get((key, norm))
            r['p_'+norm] = (math.log(old[1]/e)/math.log(int(r['N'])/old[0])
                            if old and old[1] > NOISE and e > NOISE and not constrained else '')
            previous[key, norm] = int(r['N']), e
    write(output / 'all_errors.csv', records)
    finest = [r for r in records if int(r['divisor']) == divisors[-1]]
    write(output / 'errors.csv', finest)
    write(output / 'checks.csv', checks)
    lookup = {(r['family'], int(r['N']), int(r['divisor']), r['kind'], r['group'], int(r['component'])): r for r in records}
    temporal = []
    d0, d1, d2 = divisors[-3:]
    for f in args.families:
        for n in sizes:
            fields = {d: field(Path(str(name(f, n, d))+'.field')) for d in (d0, d1, d2)}
            if any(len(a) != (n+1)**2 for a in fields.values()):
                raise ValueError('Incomplete field snapshot')
            for group in GROUPS:
                for kind in KINDS:
                    for k in range(4):
                        coarse = delta_norm(fields[d0], fields[d1], n, group, kind, k)
                        fine = delta_norm(fields[d1], fields[d2], n, group, kind, k)
                        reference = lookup[f, n, d2, kind, group, k]
                        constrained = kind.endswith('conservative') and group == 'boundary'
                        for j, norm in enumerate(NORMS):
                            error = float(reference[norm])
                            tested = kind in PRIMARY and not constrained and error > NOISE
                            ratio = fine[j]/error if error > NOISE else ''
                            temporal.append(dict(family=f, N=n, kind=kind, group=group, component=k, norm=norm,
                                divisor_coarse=d0, divisor_mid=d1, divisor_fine=d2,
                                delta_coarse=coarse[j], delta_fine=fine[j], error_fine=error,
                                fraction_error=ratio,
                                p_time=math.log2(coarse[j]/fine[j]) if min(coarse[j], fine[j]) > NOISE and not constrained else '',
                                tested=int(tested), pass_fraction=int(not tested or fine[j] <= 0.05*error)))
    write(output / 'temporal.csv', temporal)
    stability = []
    for f in args.families:
        for n in sizes[1:]:
            for kind in KINDS:
                for group in GROUPS:
                    for k in range(4):
                        for norm in NORMS:
                            a = lookup[f, n, d1, kind, group, k]['p_'+norm]
                            b = lookup[f, n, d2, kind, group, k]['p_'+norm]
                            if a == '' or b == '':
                                continue
                            change = abs(float(a)-float(b))
                            stability.append(dict(family=f, N=n, kind=kind, group=group, component=k, norm=norm,
                                divisor_mid=d1, divisor_fine=d2, p_mid=a, p_fine=b, change=change,
                                tested=int(kind in PRIMARY), pass_order=int(kind not in PRIMARY or change <= 0.1)))
    if stability:
        write(output / 'order_stability.csv', stability)
    failures = [r for r in temporal if r['tested'] and r['norm'] == 'l2' and not r['pass_fraction']]
    failures_order = [r for r in stability if r['tested'] and r['norm'] == 'l2' and not r['pass_order']]
    print(f'{len(tasks)} paired trajectories: {len(failures)} L2 temporal failures, {len(failures_order)} L2 order-stability failures', flush=True)
    for r in failures[:8]:
        print('TIME_NOT_RESOLVED', r, flush=True)
    for r in failures_order[:8]:
        print('ORDER_NOT_STABLE', r, flush=True)
    if failures or failures_order:
        raise SystemExit('Temporal independence was not demonstrated: reduce dt and repeat.')


if __name__ == '__main__':
    main()
