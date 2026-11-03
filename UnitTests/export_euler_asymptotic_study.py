#!/usr/bin/env python3
"""6B.3 signed-defect projections and fieldwise Richardson diagnostics.

Only postprocesses native harness exports. Does not extrapolate a numerical
trajectory inside either solver. Requires NumPy; output is scratch-only.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import numpy as np

GROUPS = ('all', 'boundary', 'ring1', 'ring2', 'deep3', 'fixed', 'interior')
NORMS = ('rms', 'l1', 'l2', 'linf')
KINDS = ('W_primitive', 'U_primitive', 'W_conservative', 'U_conservative', 'difference_W', 'difference_U')
FLOOR = 1e-12


def write(path, rows):
    with path.open('w') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def conservative(w):
    rho, u, v, p = w.T
    return np.column_stack((rho, rho*u, rho*v, p/0.4+0.5*rho*(u*u+v*v)))


def primitive(u):
    rho, mx, my, e = u.T
    w = np.column_stack((rho, mx/rho, my/rho, 0.4*(e-(mx*mx+my*my)/(2*rho))))
    if not np.all(np.isfinite(w)) or min(w[:, 0].min(), w[:, 3].min()) <= 0:
        raise ValueError('Diagnostic extrapolation is not admissible')
    return w


def physical(wp, uc):
    up, wc = conservative(wp), primitive(uc)
    return dict(zip(KINDS, (wp, up, wc, uc, wp-wc, up-uc)))


def masks(x, layer):
    return dict(zip(GROUPS, (np.ones(len(x), dtype=bool), layer == 0, layer == 1, layer == 2,
        layer >= 3, (x[:, 0] >= 0.25) & (x[:, 0] <= 0.75) & (x[:, 1] >= 0.25) & (x[:, 1] <= 0.75), layer >= 1)))


def norms(e, v):
    volume = v.sum()
    square = np.sum(v[:, None]*e*e, axis=0)
    return dict(rms=np.sqrt(np.mean(e*e, axis=0)), l1=np.sum(v[:, None]*np.abs(e), axis=0)/volume,
                l2=np.sqrt(square/volume), linf=np.max(np.abs(e), axis=0)), square


def rates(rows, keys, columns):
    previous = {}
    for r in sorted(rows, key=lambda x: (int(x['N']), str(x))):
        key = tuple(r[k] for k in keys)
        old = previous.get(key)
        for c in columns:
            v, ov = float(r[c]), float(old[c]) if old else 0
            r['p_'+c] = math.log(abs(ov/v))/math.log(int(r['N'])/int(old['N'])) if old and min(abs(ov), abs(v)) > FLOOR else ''
        previous[key] = r


def evolution(base, time, sizes, families, out, reference_baseline=None):
    errors, temporal, rings, checks, reproduction = [], [], [], [], []
    for family in families:
        for n in sizes:
            paths = {d: base/f'{family}_{n}_{d}.csv' for d in (1, 2, 4)}
            a = {d: np.loadtxt(str(path)+'.field') for d, path in paths.items()}
            for d in (1, 2, 4):
                if a[d].shape != ((n+1)**2, 28):
                    raise ValueError('Incomplete snapshot')
                if not np.array_equal(a[1][:, :4], a[d][:, :4]) or not np.array_equal(a[1][:, 20:], a[d][:, 20:]):
                    raise ValueError('Snapshot geometry/reference mismatch')
                with Path(str(paths[d])+'.checks').open() as stream:
                    checks.extend(dict(time=time, **r) for r in csv.DictReader(stream))
            raw = {d: physical(a[d][:, 4:8], a[d][:, 16:20]) for d in (1, 2, 4)}
            ext = physical(2*a[4][:, 4:8]-a[2][:, 4:8], 2*a[4][:, 16:20]-a[2][:, 16:20])
            old_ext = physical(2*a[2][:, 4:8]-a[1][:, 4:8], 2*a[2][:, 16:20]-a[1][:, 16:20])
            reference = {k: a[4][:, 20:24] if k.startswith('W_') else a[4][:, 24:28] if k.startswith('U_')
                         else np.zeros((len(a[4]), 4)) for k in KINDS}
            v, x = a[4][:, 1], a[4][:, 2:4]
            i, j = np.arange(len(v)) % (n+1), np.arange(len(v)) // (n+1)
            mask = masks(x, np.minimum.reduce([i, j, n-i, n-j]))
            for kind in KINDS:
                global_c = {}
                for variant, field in [('FE', raw[4]), ('Richardson', ext)]:
                    global_c[variant] = np.sum(v[:, None]*(field[kind]-reference[kind])**2, axis=0)
                    pieces = np.zeros(4)
                    for group, m in mask.items():
                        ns, c = norms(field[kind][m]-reference[kind][m], v[m])
                        if group in ('boundary', 'ring1', 'ring2', 'deep3'):
                            pieces += c
                        for k in range(4):
                            row = dict(time=time, family=family, N=n, variant=variant, kind=kind, group=group, component=k,
                                       count=int(m.sum()), volume=float(v[m].sum()), **{q: ns[q][k] for q in NORMS})
                            errors.append(row)
                            rings.append(dict(time=time, family=family, N=n, variant=variant, kind=kind, group=group,
                                component=k, volume=row['volume'], C=c[k], C_all=global_c[variant][k],
                                fraction=c[k]/global_c[variant][k] if global_c[variant][k] > 0 else 0,
                                disjoint_partition=int(group in ('boundary', 'ring1', 'ring2', 'deep3'))))
                    if not np.allclose(pieces, global_c[variant], rtol=5e-13, atol=1e-25):
                        raise ValueError('Overlapping/missing ring accounting')
                for group, m in mask.items():
                    coarse, _ = norms(raw[1][kind][m]-raw[2][kind][m], v[m])
                    fine, _ = norms(raw[2][kind][m]-raw[4][kind][m], v[m])
                    correction, _ = norms(ext[kind][m]-raw[4][kind][m], v[m])
                    uncertainty, _ = norms(ext[kind][m]-old_ext[kind][m], v[m])
                    ef, _ = norms(raw[4][kind][m]-reference[kind][m], v[m])
                    e0, _ = norms(ext[kind][m]-reference[kind][m], v[m])
                    for k in range(4):
                        for norm in NORMS:
                            d1, d2, e = coarse[norm][k], fine[norm][k], ef[norm][k]
                            constrained = kind in ('W_conservative', 'U_conservative') and group == 'boundary'
                            temporal.append(dict(time=time, family=family, N=n, kind=kind, group=group, component=k, norm=norm,
                                delta1=d1, delta2=d2, p_time=math.log2(d1/d2) if min(d1, d2) > FLOOR and not constrained else '',
                                FE_error=e, extrapolated_error=e0[norm][k], correction=correction[norm][k],
                                successive_extrapolation_difference=uncertainty[norm][k],
                                correction_fraction=correction[norm][k]/e if e > FLOOR else '',
                                extrapolation_fraction=uncertainty[norm][k]/e0[norm][k] if e0[norm][k] > FLOOR else '',
                                tested=int(kind in KINDS[:4] and not constrained and e > FLOOR)))
            if time == 0.1:
                for d in (1, 2, 4):
                    if reference_baseline is None:
                        raise ValueError('A preserved 6B.2 baseline is required for reproduction')
                    original = reference_baseline/paths[d].name
                    if not original.exists():
                        raise ValueError('Missing original baseline: '+str(original))
                    if original.exists():
                        with original.open() as stream:
                            old = list(csv.DictReader(stream))
                        with paths[d].open() as stream:
                            new = list(csv.DictReader(stream))
                        if len(old) != len(new):
                            raise ValueError('Baseline record count differs')
                        maximum = max(abs(float(b[q])-float(c[q])) for b, c in zip(old, new) for q in NORMS)
                        field_error = np.max(np.abs(np.loadtxt(str(original)+'.field')-a[d][:, :20]))
                        if max(maximum, field_error) > 2e-13:
                            raise ValueError('Baseline not reproduced')
                        reproduction.append(dict(family=family, N=n, divisor=d, max_norm_difference=maximum, max_field_difference=field_error))
    rates(errors, ('time', 'family', 'variant', 'kind', 'group', 'component'), NORMS)
    # Strong Dirichlet boundary is prescribed, not a measured spatial convergence.
    for r in errors:
        if r['kind'] in ('W_conservative', 'U_conservative') and r['group'] == 'boundary':
            for q in NORMS:
                r['p_'+q] = ''
    write(out/f'evolution_{time}.csv', errors)
    write(out/f'temporal_{time}.csv', temporal)
    write(out/f'rings_{time}.csv', rings)
    write(out/f'checks_{time}.csv', checks)
    paired = {(r['family'], r['N'], r['kind'], r['group'], r['component'], r['variant']): r for r in errors}
    orders = []
    for r in errors:
        if r['variant'] != 'FE':
            continue
        ext_row = paired[r['family'], r['N'], r['kind'], r['group'], r['component'], 'Richardson']
        for norm in NORMS:
            pf, pe = r['p_'+norm], ext_row['p_'+norm]
            orders.append(dict(time=time, family=r['family'], N=r['N'], kind=r['kind'], group=r['group'],
                component=r['component'], norm=norm, FE_error=r[norm], extrapolated_error=ext_row[norm],
                p_FE=pf, p_extrapolated=pe, order_change=float(pe)-float(pf) if pf != '' and pe != '' else ''))
    write(out/f'order_comparison_{time}.csv', orders)
    if reproduction:
        write(out/'reproduction.csv', reproduction)


def defects(base, out):
    records, moments, parts, checks = [], [], [], []
    raw_export = out/'signed_defects.csv.gz'
    with gzip.open(raw_export, 'wt') as stream:
        writer = csv.writer(stream)
        writer.writerow(('family', 'N', 'time', 'node', 'x', 'y', 'volume', 'layer', 'free', 'component', 'tau_W', 'J_tau_W', 'tau_U_free'))
        for family in ('regular', 'smooth', 'jitter'):
            for n in (8, 16, 32, 64, 128):
                source = base/f'{family}_{n}_0.csv'
                with source.open() as f:
                    checks.extend(r for r in csv.DictReader(f) if r['group'] == 'check')
                for t in (0., .05, .1):
                    a = np.loadtxt(str(source)+f'.signed_{t:.6f}')
                    if a.shape != ((n+1)**2, 37) or not np.all(np.isfinite(a)):
                        raise ValueError('Incomplete signed export')
                    x, v, layer = a[:, 1:3], a[:, 3], a[:, 4]
                    mask = masks(x, layer)
                    for row in a:
                        for k in range(4):
                            writer.writerow((family, n, t, int(row[0]), row[1], row[2], row[3], int(row[4]), int(row[4] > 0), k,
                                             row[5+k], row[9+k], row[13+k] if row[4] > 0 else ''))
                    for kind, data in [('tau_W', a[:, 5:9]), ('J_tau_W', a[:, 9:13]), ('tau_U', a[:, 13:17]),
                                       ('difference_U', a[:, 9:13]-a[:, 13:17])]:
                        free = kind != 'tau_W'
                        total = np.sum(v[(layer > 0) if free else np.ones(len(v), dtype=bool), None]*
                            data[(layer > 0) if free else np.ones(len(v), dtype=bool)]**2, axis=0)
                        partition = np.zeros(4)
                        for group, original_mask in mask.items():
                            m = original_mask & (layer > 0) if free else original_mask
                            if not m.any():
                                continue
                            ns, c = norms(data[m], v[m])
                            if group in ('boundary', 'ring1', 'ring2', 'deep3'):
                                partition += c
                            for k in range(4):
                                records.append(dict(time=t, family=family, N=n, kind=kind, group='all_free' if free and group == 'all' else group,
                                    component=k, count=int(m.sum()), volume=v[m].sum(), C=c[k], C_total=total[k],
                                    fraction=c[k]/total[k] if total[k] else 0, **{q: ns[q][k] for q in NORMS}))
                        if not np.allclose(partition, total, rtol=5e-13, atol=1e-25):
                            raise ValueError('Defect ring accounting failed')
                        m = mask['fixed']
                        xm, ym = x[m].T
                        phi = {'1': np.ones(m.sum()), 'x': xm, 'y': ym,
                            'sin_sin': np.sin(2*np.pi*xm)*np.sin(2*np.pi*ym),
                            'cos_cos': np.cos(2*np.pi*xm)*np.cos(2*np.pi*ym)}
                        for label, values in phi.items():
                            signed = np.sum(v[m, None]*values[:, None]*data[m], axis=0)/0.25
                            absolute_integral = np.sum(v[m, None]*np.abs(values[:, None]*data[m]), axis=0)/0.25
                            for k in range(4):
                                moments.append(dict(time=t, family=family, N=n, kind=kind, group='fixed', component=k,
                                    phi=label, normalization_area=0.25, sampled_volume=v[m].sum(), signed=signed[k],
                                    absolute=abs(signed[k]), integral_absolute=absolute_integral[k],
                                    cancellation_ratio=abs(signed[k])/absolute_integral[k] if absolute_integral[k] else 0))
                    # Native API does not expose the same fluctuation/volume split.
                    for kind, start in [('faces', 17), ('volume', 21), ('boundary', 25), ('source', 29), ('rhs', 33)]:
                        data = a[:, start:start+4]
                        for group, m in mask.items():
                            ns, c = norms(data[m], v[m])
                            for k in range(4):
                                parts.append(dict(time=t, family=family, N=n, kind=kind, group=group, component=k,
                                    volume=v[m].sum(), C=c[k], **{q: ns[q][k] for q in NORMS}))
    rates(records, ('time', 'family', 'kind', 'group', 'component'), NORMS)
    rates(moments, ('time', 'family', 'kind', 'group', 'component', 'phi'), ('absolute', 'integral_absolute'))
    write(out/'defects.csv', records)
    write(out/'moments.csv', moments)
    write(out/'contributions.csv', parts)
    write(out/'defect_checks.csv', checks)
    (out/'signed_manifest.json').write_text(json.dumps(dict(file=raw_export.name,
        sha256=hashlib.sha256(raw_export.read_bytes()).hexdigest(), representation='pointwise primitive / matched conservative',
        conservative_boundary='excluded; no zero truncation error inferred'), indent=2)+'\n')


def package_outputs(root, out, preservation):
    for kind in ('evolution', 'temporal', 'rings'):
        rows = []
        for t in ('0.1', '0.2'):
            with (out/f'{kind}_{t}.csv').open() as stream:
                rows.extend(csv.DictReader(stream))
        write(out/f'{kind}.csv', rows)
    rows = []
    for t in ('0.1', '0.2'):
        with (out/f'checks_{t}.csv').open() as stream:
            rows.extend(dict(stage='evolution', **r) for r in csv.DictReader(stream))
    with (out/'defect_checks.csv').open() as stream:
        for r in csv.DictReader(stream):
            rows.append(dict(stage='defect', time=r['time'], family=r['family'], N=r['N'], divisor='NA',
                             kind=r['kind'], component=r['component'], value=r['rms']))
    write(out/'checks.csv', rows)
    manifest = json.loads((out/'signed_manifest.json').read_text())
    manifest['raw_location'] = str((out/'signed_defects.csv.gz').resolve())
    manifest['prior_files_sha256'] = json.loads(preservation.read_text()) if preservation else {}
    for path, expected in manifest['prior_files_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError('Prior file changed: '+path)
    manifest['builds'] = {k: json.loads((root/k/'build_manifest.json').read_text())
                          for k in ('defect-build', 'evolution-build', 'long-build')}
    (out/'provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--skip-defects', action='store_true', help='Use only after a completed signed-defect export')
    p.add_argument('--reference-baseline', type=Path, default=Path('/private/tmp/su2-mms-6b2-complete/campaign'))
    p.add_argument('--preservation-manifest', type=Path, help='Pre-run mapping of prior file paths to SHA256')
    a = p.parse_args()
    out = a.output_dir
    out.mkdir(parents=True, exist_ok=True)
    evolution(a.root/'baseline', .1, [8, 16, 32, 64], ['regular', 'smooth', 'jitter'], out, a.reference_baseline)
    evolution(a.root/'long', .2, [8, 16, 32], ['regular', 'jitter'], out)
    if not a.skip_defects:
        defects(a.root/'defects', out)
    package_outputs(a.root, out, a.preservation_manifest)
    print('Evolution, temporal extrapolation, signed defects, moments and disjoint ledgers exported.')


if __name__ == '__main__':
    main()
