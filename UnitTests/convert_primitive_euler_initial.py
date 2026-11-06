#!/usr/bin/env python3
"""Convert a same-mesh SU2 Euler Float64 restart to primitive FILE initialization.

Binary SU2 stores rows in implicit GlobalIndex order. This is verified against
the original mesh coordinates; no interpolation or nearest-neighbor matching.
The input restart is never written. Output is exclusive-create.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import numpy as np
from run_primitive_euler_driver_e2e import binary


def mesh_coordinates(mesh):
    lines=mesh.read_text().splitlines()
    at=next(i for i,line in enumerate(lines) if line.startswith('NPOIN='))
    n=int(lines[at].split('=')[1].split()[0]); xy=np.empty((n,2)); seen=set()
    for line in lines[at+1:at+1+n]:
        fields=line.split(); index=int(fields[2])
        if index in seen or not 0<=index<n: raise ValueError('Duplicate/out-of-range mesh GlobalIndex')
        seen.add(index); xy[index]=[float(x) for x in fields[:2]]
    if len(seen)!=n or not np.all(np.isfinite(xy)): raise ValueError('Incomplete/nonfinite mesh coordinates')
    return xy


def write_initial(path,xy,w):
    names=['x','y','Density','Velocity_x','Velocity_y','Pressure','GlobalIndex']
    n=len(xy)
    if w.shape!=(n,4) or not np.all(np.isfinite(w)) or np.min(w[:,[0,3]])<=0:
        raise ValueError('Invalid primitive initial field')
    values=np.column_stack([xy,w,np.arange(n)]).astype(np.float64)
    with path.open('xb') as f:
        f.write(struct.pack('=5i',535532,7,n,0,0))
        for name in names: f.write(name.encode().ljust(33,b'\0'))
        f.write(values.tobytes())


def convert(source,mesh,target,gamma=1.4):
    if not np.isfinite(gamma) or gamma<=1: raise ValueError('Invalid gamma')
    checksum=hashlib.sha256(source.read_bytes()).hexdigest()
    data=binary(source); xy=mesh_coordinates(mesh); n=len(xy)
    expected=['x','y','Density','Momentum_x','Momentum_y','Energy']
    if list(data)!=expected or any(len(data[k])!=n for k in expected):
        raise ValueError('Expected complete 2D Euler conservative binary restart')
    coord=np.column_stack([data['x'],data['y']])
    np.testing.assert_allclose(coord,xy,rtol=0,atol=1e-12)
    rho=data['Density']; momentum=np.column_stack([data['Momentum_x'],data['Momentum_y']])
    if np.min(rho)<=0: raise ValueError('Nonpositive native density')
    velocity=momentum/rho[:,None]
    pressure=(gamma-1)*(data['Energy']-.5*np.sum(momentum**2,axis=1)/rho)
    w=np.column_stack([rho,velocity,pressure])
    inverse=np.column_stack([rho,rho[:,None]*velocity,pressure/(gamma-1)+.5*rho*np.sum(velocity**2,axis=1)])
    original=np.column_stack([rho,momentum,data['Energy']])
    np.testing.assert_allclose(inverse,original,rtol=8*np.finfo(float).eps,atol=0)
    write_initial(target,xy,w)
    reloaded=binary(target)
    np.testing.assert_array_equal(reloaded['GlobalIndex'],np.arange(n))
    np.testing.assert_array_equal(np.column_stack([reloaded[k] for k in ['Density','Velocity_x','Velocity_y','Pressure']]),w)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==checksum
    return dict(source=str(source),source_sha256=checksum,mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),
                nodes=n,gamma=gamma,coordinate_max_error=float(np.max(abs(coord-xy))),
                inverse_max_error=float(np.max(abs(inverse-original))),min_rho=float(np.min(rho)),
                min_p=float(np.min(pressure)),identity='implicit binary GlobalIndex verified against explicit mesh IDs',
                projection='none in converter; explicit solver initialization projection is measured separately')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True); p.add_argument('--mesh',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); p.add_argument('--gamma',type=float,default=1.4)
    a=p.parse_args(); print(json.dumps(convert(a.source,a.mesh,a.output,a.gamma),indent=2))
