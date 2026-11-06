#!/usr/bin/env python3
"""Bounded phase-15 NACA AOA=0 campaign using the actual SU2_CFD drivers/writers.

No numerical solver or VTK writer is implemented here. Native initial/final
audits use single explicit iterations at CFL=0; native wall reprojection is
measured separately from the zero temporal increment, on the same mesh.
Use --postprocess-only to regenerate plots/checks without rerunning either solver.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import time
import xml.etree.ElementTree as ET

import numpy as np
from run_primitive_euler_driver_e2e import binary

ROOT = Path(__file__).resolve().parents[1]
CASES = ROOT / 'TestCases/primitive_euler_examples/naca0012_aoa0'
MESH = ROOT / 'QuickStart/mesh_NACA0012_inv.su2'


def table(path):
    with path.open() as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        reader.fieldnames = [k.strip() for k in reader.fieldnames]
        return [{k: float(v) for k, v in row.items()} for row in reader]


def write(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def setting(text, key, value):
    line = f'{key}={value}'
    return re.sub(rf'(?m)^{key}\s*=.*$', line, text) if re.search(rf'(?m)^{key}\s*=', text) else text + line + '\n'


def configuration(solver, steps):
    text = (CASES / f'{solver}.cfg').read_text()
    text = setting(text, 'MESH_FILENAME', str(MESH))
    text = setting(text, 'TIME_ITER' if solver == 'primitive' else 'ITER', steps)
    frequency = max(1, steps // 5)
    for key in ['OUTPUT_WRT_FREQ', 'SCREEN_WRT_FREQ_TIME' if solver == 'primitive' else 'SCREEN_WRT_FREQ_INNER']:
        text = setting(text, key, frequency)
    return text


def run(exe, folder, text, seconds, expected=0):
    folder.mkdir(parents=True, exist_ok=False)
    (folder / 'case.cfg').write_text(text)
    start = time.monotonic()
    timeout = False
    with (folder / 'run.log').open('w') as log:
        process = subprocess.Popen([str(exe), 'case.cfg'], cwd=folder, stdout=log,
                                   stderr=subprocess.STDOUT, env={**os.environ, 'OMP_NUM_THREADS': '1'})
        try:
            process.wait(timeout=seconds)
        except subprocess.TimeoutExpired:
            timeout = True
            process.terminate()  # SU2 handles SIGTERM and writes the last accepted state.
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
    result = dict(seconds=time.monotonic()-start, timeout=timeout, returncode=process.returncode)
    (folder / 'execution.json').write_text(json.dumps(result, indent=2))
    if expected == 0 and process.returncode != 0:
        raise RuntimeError(f'{folder}: solver failed; inspect run.log')
    if expected != 0 and process.returncode == 0:
        raise AssertionError(f'{folder}: invalid input was accepted')
    return result


def vtu(path):
    """Validate/read SU2's raw-appended little-endian VTU, including byte counts."""
    raw = path.read_bytes()
    begin = raw.index(b'<AppendedData')
    header = ET.fromstring(raw[:begin] + b'</VTKFile>')
    assert header.attrib['byte_order'] == 'LittleEndian' and header.attrib['header_type'] == 'UInt64'
    data_begin = raw.index(b'_', begin) + 1
    assert raw.rstrip().endswith(b'</VTKFile>')
    piece = header.find('UnstructuredGrid/Piece')
    points, cells = int(piece.attrib['NumberOfPoints']), int(piece.attrib['NumberOfCells'])
    types = {'Float32': '<f4', 'Float64': '<f8', 'Int32': '<i4', 'UInt8': 'u1'}

    def array(node):
        at = data_begin + int(node.attrib['offset'])
        size, = struct.unpack_from('<Q', raw, at)
        dtype = np.dtype(types[node.attrib['type']])
        assert at + 8 + size <= len(raw) and size % dtype.itemsize == 0
        values = np.frombuffer(raw, dtype=dtype, count=size // dtype.itemsize, offset=at + 8).copy()
        assert np.all(np.isfinite(values))
        return values.reshape(-1, int(node.attrib['NumberOfComponents'])).squeeze(axis=1) if int(node.attrib['NumberOfComponents']) == 1 else values.reshape(-1, int(node.attrib['NumberOfComponents']))

    coordinates = array(piece.find('Points/DataArray'))
    assert coordinates.shape == (points, 3)
    topology = {node.attrib['Name']: array(node) for node in piece.findall('Cells/DataArray')}
    assert len(topology['offsets']) == cells and len(topology['types']) == cells
    assert topology['offsets'][-1] == len(topology['connectivity'])
    assert np.min(topology['connectivity']) >= 0 and np.max(topology['connectivity']) < points
    fields = {node.attrib['Name']: array(node) for node in piece.findall('PointData/DataArray')}
    assert all(len(values) == points for values in fields.values())
    return coordinates, fields, cells


def snapshots(folder, solver):
    records = []
    mesh = MESH.read_text()
    expected_points = int(re.search(r'NPOIN=\s*(\d+)', mesh)[1])
    expected_cells = int(re.search(r'NELEM=\s*(\d+)', mesh)[1])
    for path in sorted(folder.glob(f'{solver}_*.vtm')):
        root = ET.parse(path).getroot()
        datasets = root.findall('.//DataSet')
        for node in datasets:
            target = path.parent / node.attrib['file']
            assert target.is_file() and target.stat().st_size > 0
            xyz, fields, cells = vtu(target)
            if node.attrib['name'] != 'Internal':
                continue
            # Some native steady VTM files contain both unindexed and indexed
            # blocks. Validate all files, but never mix different snapshots.
            if solver == 'native' and Path(node.attrib['file']).parts[0] != path.stem:
                continue
            assert len(xyz) == expected_points and cells == expected_cells
            rho, p, velocity = fields['Density'], fields['Pressure'], fields['Velocity']
            if solver == 'native':
                # Check that displayed primitives and current conservative fields agree,
                # allowing only native Float32 output precision.
                momentum = fields['Momentum'][:, :2].astype(float)
                density = rho.astype(float)
                physical_p = .4*(fields['Energy'].astype(float)-.5*np.sum(momentum*momentum, axis=1)/density)
                np.testing.assert_allclose(p, physical_p, rtol=1e-6, atol=1e-7)
                np.testing.assert_allclose(velocity[:, :2], momentum/density[:, None], rtol=1e-6, atol=1e-7)
            mach = np.linalg.norm(velocity[:, :2], axis=1) / np.sqrt(1.4*p/rho)
            assert np.min(rho) > 0 and np.min(p) > 0 and np.max(mach) < 1
            np.testing.assert_allclose(fields['Mach'], mach, rtol=5e-7, atol=5e-8)
            step = int(path.stem.rsplit('_', 1)[1])
            records.append(dict(iteration=step, file=path.name, sequence_file=node.attrib['file'], points=len(xyz), cells=cells,
                                min_rho=float(np.min(rho)), min_p=float(np.min(p)), max_mach=float(np.max(mach))))
    assert records, f'No validated snapshots in {folder}'
    return records


def pvd(folder, records, initial=None):
    root = ET.Element('VTKFile', type='Collection', version='0.1', byte_order='LittleEndian')
    collection = ET.SubElement(root, 'Collection')
    if initial:
        ET.SubElement(collection, 'DataSet', timestep='0', group='', part='0', file=initial)
        assert (folder / initial).is_file()
    for row in records:
        # Index is the abscissa for steady comparison, NEVER physical pseudo-time.
        index = row['iteration'] + (1 if initial else 0)
        target = row['sequence_file']
        ET.SubElement(collection, 'DataSet', timestep=str(index), group='', part='0', file=target)
        assert (folder / target).is_file()
    manifest = folder / 'iterations.pvd'
    ET.ElementTree(root).write(manifest, encoding='utf-8', xml_declaration=True)
    entries = ET.parse(manifest).getroot().findall('Collection/DataSet')
    assert np.all(np.diff([float(node.attrib['timestep']) for node in entries]) > 0)
    assert all((folder / node.attrib['file']).is_file() for node in entries)


def cp_curves(output, primitive_step, native_file):
    with (output / 'primitive' / f'surface_primitive_pressure_{primitive_step}.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    primitive = np.array([[float(r[k]) for k in ['x', 'y', 'Cp']] for r in rows])
    root = ET.parse(native_file).getroot()
    boundary = next(n for n in root.findall('.//DataSet') if n.attrib['name'] == 'airfoil'
                    and Path(n.attrib['file']).parts[0] == native_file.stem)
    xyz, fields, _ = vtu(native_file.parent / boundary.attrib['file'])
    native = np.column_stack([xyz[:, :2], (fields['Pressure'].astype(float)-1)/.028])
    differences, curves = [], []
    for side, sign in [('upper', 1), ('lower', -1)]:
        a = primitive[sign*primitive[:, 1] > 0]
        b = native[sign*native[:, 1] >= 0]
        a, b = a[np.argsort(a[:, 0])], b[np.argsort(b[:, 0])]
        assert np.all(np.diff(b[:, 0]) > 0), 'Ambiguous native Cp interpolation abscissae'
        keep = (a[:, 0] >= b[0, 0]) & (a[:, 0] <= b[-1, 0])
        a = a[keep]
        reference = np.interp(a[:, 0], b[:, 0], b[:, 2])
        diff = a[:, 2]-reference
        differences.extend(diff)
        curves.append((side, a, b))
    return curves, dict(points=len(differences), rms=float(np.sqrt(np.mean(np.square(differences)))),
                        maximum=float(np.max(np.abs(differences))))


def postprocess(out):
    assert hashlib.sha256(MESH.read_bytes()).hexdigest() == (out / 'mesh.sha256').read_text().strip()
    os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/primitive-phase15-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    compare = out / 'comparison'
    compare.mkdir(exist_ok=True)
    histories = {s: table(out / s / 'history.csv') for s in ['primitive', 'native']}
    for row in histories['primitive']:
        assert row['Min_rho'] > 0 and row['Min_p'] > 0
        assert 0 <= row['MAX_MACH'] < 1 and 0 <= row['MAX_TRACE_MACH'] < 1
    checked = {s: snapshots(out / s, s) for s in histories}
    initial = snapshots(out / 'native_initial', 'native')
    assert len(checked['primitive']) >= 4 and len(checked['native']) >= 3
    pvd(out / 'primitive', checked['primitive'])
    pvd(out / 'native', checked['native'], '../native_initial/'+initial[0]['sequence_file'])
    summary = []
    native_final = table(out / 'native_recheck' / 'history.csv')[-1]
    # CFL=0 suppresses the update, but native BCs can reproject wall momentum.
    native_checkpoint = sorted((out / 'native').glob('native_state*.dat'))[-1]
    audit_checkpoint = sorted((out / 'native_recheck').glob('native_state*.dat'))[-1]
    saved, audited = binary(native_checkpoint), binary(audit_checkpoint)
    assert saved.keys() == audited.keys()
    audit_change = max(float(np.max(np.abs(saved[k]-audited[k]))) for k in saved)
    mesh_lines = MESH.read_text().splitlines()
    point_start = next(i for i, line in enumerate(mesh_lines) if line.startswith('NPOIN='))
    count = int(mesh_lines[point_start].split('=')[1].split()[0])
    positions = [tuple(map(float, line.split()[:2])) for line in mesh_lines[point_start+1:point_start+1+count]]
    airfoil_start = next(i for i, line in enumerate(mesh_lines) if line.strip() == 'MARKER_TAG= airfoil')
    edges = int(mesh_lines[airfoil_start+1].split('=')[1])
    wall_positions = {positions[int(index)] for line in mesh_lines[airfoil_start+2:airfoil_start+2+edges]
                      for index in line.split()[1:3]}
    wall = np.array([(x, y) in wall_positions for x, y in zip(saved['x'], saved['y'])])
    for key in saved:
        difference = np.abs(saved[key]-audited[key])
        if key in ['Momentum_x', 'Momentum_y']:
            assert np.max(difference[~wall]) == 0
        else:
            assert np.max(difference) == 0
    assert audit_change < 1e-9, f'Unexpectedly large native wall reprojection: {audit_change}'
    (compare / 'native_recheck.json').write_text(json.dumps(dict(max_wall_momentum_reprojection=audit_change,
        density_energy_interior_unchanged=True,
        residual=10**native_final['rms[Rho]']), indent=2))
    for solver, history in histories.items():
        first, last = history[0], history[-1]
        residual_initial = first['RHS_CURRENT'] if solver == 'primitive' else 10**first['rms[Rho]']
        residual_final = last['RHS_CURRENT'] if solver == 'primitive' else 10**native_final['rms[Rho]']
        execution = json.loads((out / solver / 'execution.json').read_text())
        window = history[-max(2, min(100, len(history)//10)):]
        coefficients = last if solver == 'primitive' else native_final
        record = dict(solver=solver, iterations=int(last['Time_Iter'] if solver == 'primitive' else last['Inner_Iter'])+(solver == 'native'),
                      residual_initial=residual_initial, residual_final=residual_final,
                      residual_ratio=residual_final/residual_initial, reduction_decades=np.log10(residual_initial/residual_final),
                      CL=coefficients['CL'], CD=coefficients['CD'], max_mach=max(r['max_mach'] for r in checked[solver]),
                      min_rho=min(r['min_rho'] for r in checked[solver]), min_p=min(r['min_p'] for r in checked[solver]),
                      seconds=execution['seconds'], budget_stop=execution['timeout'], final_residual_recomputed=True,
                      CL_window_range=max(r['CL'] for r in window)-min(r['CL'] for r in window),
                      CD_window_range=max(r['CD'] for r in window)-min(r['CD'] for r in window))
        if solver == 'primitive':
            record['max_mach'] = max(r['MAX_MACH'] for r in history)
            record['min_rho'] = min(r['Min_rho'] for r in history)
            record['min_p'] = min(r['Min_p'] for r in history)
        summary.append(record)
        write(compare / f'{solver}_snapshots.csv', checked[solver])
    write(compare / 'summary.csv', summary)
    te_rows = []
    for path in sorted((out / 'primitive').glob('surface_primitive_te_*.csv'),
                       key=lambda p: int(p.stem.rsplit('_', 1)[1])):
        faces = table(path)
        assert len(faces) == 2 and faces[0]['global_id'] == faces[1]['global_id'] == 199
        assert (faces[0]['x'], faces[0]['y']) != (faces[1]['x'], faces[1]['y'])
        assert faces[0]['rho'] > 0 and faces[0]['p'] > 0
        te_rows.append(dict(iteration=faces[0]['iteration'], time=faces[0]['time'],
            rho=faces[0]['rho'], u=faces[0]['u'], v=faces[0]['v'], p=faces[0]['p'],
            rhs_norm=np.linalg.norm([faces[0][k] for k in ['rhs_rho', 'rhs_u', 'rhs_v', 'rhs_p']]),
            trace_un_1=faces[0]['trace_un'], trace_un_2=faces[1]['trace_un'],
            minus_norm_1=np.linalg.norm([faces[0][k] for k in ['minus_rho', 'minus_u', 'minus_v', 'minus_p']]),
            minus_norm_2=np.linalg.norm([faces[1][k] for k in ['minus_rho', 'minus_u', 'minus_v', 'minus_p']])))
    assert len(te_rows) == len(checked['primitive'])
    write(compare / 'te_summary.csv', te_rows)
    for metric, filename, ylabel in [('residual', 'residuals.png', 'Residual / initial residual'),
                                     ('CL', 'lift.png', 'CL'), ('CD', 'drag.png', 'CD')]:
        fig, ax = plt.subplots()
        for solver, history in histories.items():
            x = [r['Time_Iter'] if solver == 'primitive' else r['Inner_Iter']+1 for r in history]
            if metric == 'residual':
                y = np.array([r['RHS_CURRENT'] if solver == 'primitive' else 10**r['rms[Rho]'] for r in history])
                y /= y[0]
                ax.semilogy(x, y, label=solver)
            else:
                ax.plot(x, [r[metric] for r in history], label=solver)
        ax.set(xlabel='Accepted global FE steps / native steady iterations', ylabel=ylabel)
        ax.legend(); ax.grid(); fig.tight_layout(); fig.savefig(compare / filename); plt.close(fig)
    curves, difference = cp_curves(out, checked['primitive'][-1]['iteration'], out / 'native' / checked['native'][-1]['file'])
    write(compare / 'cp_difference.csv', [difference])
    fig, ax = plt.subplots()
    cp_rows = []
    for side, a, b in curves:
        ax.plot(a[:, 0], a[:, 2], label=f'primitive {side}')
        ax.plot(b[:, 0], b[:, 2], '--', label=f'native {side}')
        for x, y, cp in a:
            cp_rows.append(dict(side=side, x=x, y=y, primitive_cp=cp,
                                native_cp_interpolated=np.interp(x, b[:, 0], b[:, 2])))
    write(compare / 'cp_comparison.csv', cp_rows)
    ax.invert_yaxis(); ax.set(xlabel='x/c', ylabel='Cp'); ax.legend(); ax.grid()
    fig.tight_layout(); fig.savefig(compare / 'cp.png'); plt.close(fig)
    print(json.dumps(dict(summary=summary, cp=difference), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, default=ROOT / 'build/SU2_CFD/src/SU2_CFD')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--primitive-steps', type=int, default=3000)
    parser.add_argument('--native-steps', type=int, default=600)
    parser.add_argument('--seconds-per-run', type=float, default=180)
    parser.add_argument('--postprocess-only', action='store_true')
    args = parser.parse_args()
    out, exe = args.output_dir.resolve(), args.executable.resolve()
    if not args.postprocess_only:
        if args.primitive_steps < 10 or args.native_steps < 10 or args.seconds_per_run <= 0:
            parser.error('Use at least ten steps/iterations and a positive execution budget')
        out.mkdir(parents=True, exist_ok=False)
        checksum = hashlib.sha256(MESH.read_bytes()).hexdigest()
        for solver in ['primitive', 'native']:
            run(exe, out / f'smoke_{solver}', configuration(solver, 5), args.seconds_per_run)
            snapshots(out / f'smoke_{solver}', solver)
        # The explicitly selected TE must keep the freestream velocity at initialization.
        initial_vtm = sorted((out / 'smoke_primitive').glob('primitive_*.vtm'))[0]
        internal = next(n for n in ET.parse(initial_vtm).getroot().findall('.//DataSet') if n.attrib['name'] == 'Internal')
        xyz, fields, _ = vtu(initial_vtm.parent / internal.attrib['file'])
        tip = np.flatnonzero((xyz[:, 0] == 1) & (xyz[:, 1] == 0))
        assert len(tip) == 1
        np.testing.assert_allclose(fields['Velocity'][tip[0], :2], [.2*np.sqrt(1.4), 0], rtol=3e-7, atol=1e-12)
        text = configuration('primitive', 5)
        invalid = re.sub(r'(?m)^PRIMITIVE_WEAK_TE_GLOBAL_ID=.*\n', '', text)
        run(exe, out / 'reject_without_optin', invalid, args.seconds_per_run, expected=1)
        assert 'trailing edge' in (out / 'reject_without_optin/run.log').read_text().lower()
        run(exe, out / 'reject_missing_id', setting(text, 'PRIMITIVE_WEAK_TE_GLOBAL_ID', 999999999), args.seconds_per_run, expected=1)
        run(exe, out / 'reject_smooth_wall', setting(text, 'PRIMITIVE_WEAK_TE_GLOBAL_ID', 50), args.seconds_per_run, expected=1)
        run(exe, out / 'reject_native_option', setting(configuration('native', 5), 'PRIMITIVE_WEAK_TE_GLOBAL_ID', 199), args.seconds_per_run, expected=1)
        initial = configuration('native', 1)
        for key, value in [('CFL_NUMBER', 0), ('MGLEVEL', 0), ('TIME_DISCRE_FLOW', 'EULER_EXPLICIT')]:
            initial = setting(initial, key, value)
        run(exe, out / 'native_initial', initial, args.seconds_per_run)
        run(exe, out / 'primitive', configuration('primitive', args.primitive_steps), args.seconds_per_run)
        run(exe, out / 'native', configuration('native', args.native_steps), args.seconds_per_run)
        # Re-evaluate native final spatial residual without advancing/replacing its scheme.
        checkpoints = sorted((out / 'native').glob('native_state*.dat'))
        assert checkpoints
        recheck = setting(initial, 'RESTART_SOL', 'YES')
        recheck = setting(recheck, 'SOLUTION_FILENAME', str(checkpoints[-1]))
        run(exe, out / 'native_recheck', recheck, args.seconds_per_run)
        assert hashlib.sha256(MESH.read_bytes()).hexdigest() == checksum
        (out / 'mesh.sha256').write_text(checksum+'\n')
    postprocess(out)


if __name__ == '__main__':
    main()
