#!/usr/bin/env python3
"""Small polygonal-wall driver smoke; no convergence campaign or native oracle."""
import argparse
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from primitive_euler_driver_cases import config
from run_primitive_euler_driver_e2e import binary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    xy = np.array([[0, 0], [1, .2], [2, .3], [2, 1.3], [0, 1], [.9, .7]])
    cells = [(0, 1, 5), (1, 2, 5), (2, 3, 5), (3, 4, 5), (4, 0, 5)]
    mesh = ['NDIME=2', 'NELEM=5'] + ['5 ' + ' '.join(map(str, t)) for t in cells]
    mesh += ['NPOIN=6'] + [f'{x:.17g} {y:.17g} {i}' for i, (x, y) in enumerate(xy)]
    mesh += ['NMARK=2', 'MARKER_TAG=wall', 'MARKER_ELEMS=2', '3 0 1', '3 1 2',
             'MARKER_TAG=open', 'MARKER_ELEMS=3', '3 2 3', '3 3 4', '3 4 0']
    (out / 'polygon.su2').write_text('\n'.join(mesh) + '\n')
    cfg = config('polygon.su2', 8, uniform=True, steps=5, dt=1e-4)
    cfg = cfg.replace('PRIMITIVE_INITIALIZATION=UNIFORM', 'PRIMITIVE_INITIALIZATION=FREESTREAM')
    cfg = cfg.replace('MARKER_CUSTOM=(boundary)', 'MARKER_EULER=(wall)\nMARKER_CUSTOM=(open)')
    cfg = cfg.replace('PRIMITIVE_BOUNDARY_STATE=(2,0.5,0.3,3)', 'PRIMITIVE_BOUNDARY_STATE=(1,0.4,0.2,1)')
    cfg = cfg.replace('MACH_NUMBER=0.4', f'MACH_NUMBER={math.sqrt(.2/1.4):.17g}')
    cfg += f'GAS_CONSTANT=1\nAOA={math.degrees(math.atan2(.2,.4)):.17g}\n'
    (out / 'case.cfg').write_text(cfg)
    result = subprocess.run([str(args.executable.resolve()), 'case.cfg'], cwd=out,
                            env={**os.environ, 'OMP_NUM_THREADS': '1'},
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (out / 'run.log').write_text(result.stdout)
    if result.returncode or 'tangential velocity projection at 3 Euler-wall nodes' not in result.stdout:
        raise RuntimeError(str(out / 'run.log'))
    normals = np.array([[.1, -.5], [.15, -1], [.05, -.5]])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    snapshots = sorted(out.glob('primitive_state_*.dat'))
    if len(snapshots) < 2:
        raise AssertionError('Expected initial and evolved restart snapshots')
    for path in snapshots:
        data = binary(path)
        # Locate by coordinates: the production reader may reorder mesh nodes.
        coordinates = np.column_stack([data['x'], data['y']])
        indices = [int(np.argmin(np.linalg.norm(coordinates - p, axis=1))) for p in xy]
        np.testing.assert_allclose(coordinates[indices], xy, atol=1e-13, rtol=0)
        velocity = np.column_stack([data['Velocity_x'], data['Velocity_y']])[indices]
        assert np.max(np.abs(np.sum(velocity[:3] * normals, axis=1))) < 2e-13
        assert np.min(data['Density']) > 0 and np.min(data['Pressure']) > 0
        if path == snapshots[0]:
            np.testing.assert_allclose(data['Density'], 1, atol=1e-14, rtol=0)
            np.testing.assert_allclose(data['Pressure'], 1, atol=1e-14, rtol=0)
            np.testing.assert_allclose(velocity[3:], np.tile([.4, .2], (3, 1)), atol=1e-14, rtol=0)
            assert np.min(np.linalg.norm(velocity[:3], axis=1)) > .3
    print(f'PASS polygon CLI: {len(snapshots)} finite positive snapshots, nodal slip preserved')


if __name__ == '__main__':
    main()
