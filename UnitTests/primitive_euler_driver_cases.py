#!/usr/bin/env python3
"""Generate ordinary SU2 input files from the exact archived phase-7 mesh.

Only stdlib is needed to generate examples. No solver harness is executed.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]


def mesh(family,n):
    if family not in ('regular','jitter') or n not in (8,16,32,64):
        raise ValueError('Use regular/jitter and N=8,16,32,64')
    path=ROOT/f'Docs/primitive_euler_phase7_vortex/{family}_{n}_0.vtu'
    piece=ET.parse(path).getroot().find('UnstructuredGrid/Piece')
    values=list(map(float,piece.find('Points/DataArray').text.split()))
    xy=[values[i:i+2] for i in range(0,len(values),3)]
    cells={a.attrib['Name']:list(map(int,a.text.split())) for a in piece.findall('Cells/DataArray')}
    if set(cells['types'])!={5}:raise ValueError('Expected triangular archived mesh')
    triangles=[cells['connectivity'][i:i+3] for i in range(0,len(cells['connectivity']),3)]
    if len(xy)!=(n+1)**2 or len(triangles)!=2*n*n:raise ValueError('Bad archived dimensions')
    edges=defaultdict(list)
    for t in triangles:
        for k in range(3):edges[tuple(sorted((t[k],t[(k+1)%3])))].append((t[k],t[(k+1)%3]))
    boundary=[pairs[0] for _,pairs in sorted(edges.items()) if len(pairs)==1]
    text=['NDIME= 2',f'NELEM= {len(triangles)}']
    text.extend('5 '+' '.join(map(str,t)) for t in triangles)
    text.append(f'NPOIN= {len(xy)}')
    text.extend(f'{x:.17g} {y:.17g} {i}' for i,(x,y) in enumerate(xy))
    text.extend(['NMARK= 1','MARKER_TAG= boundary',f'MARKER_ELEMS= {len(boundary)}'])
    text.extend(f'3 {i} {j}' for i,j in boundary)
    return '\n'.join(text)+'\n'


def config(mesh_name,n,uniform=False,steps=None,dt=None,prefix=None):
    count=steps or (8 if uniform else 2*n*n)  # phase-7 finest divisor d=8
    step=dt if dt is not None else (.01 if uniform else .5/count)
    mode='UNIFORM' if uniform else 'ISENTROPIC_VORTEX'
    text=f'''% Native single-zone primitive solver; no source, no general conservative BC.
SOLVER=PRIMITIVE_EULER
MATH_PROBLEM=DIRECT
MESH_FILENAME={mesh_name}
MESH_FORMAT=SU2
TIME_DOMAIN=YES
TIME_MARCHING=TIME_STEPPING
TIME_DISCRE_FLOW=EULER_EXPLICIT
TIME_STEP={step:.17g}
TIME_ITER={count}
MAX_TIME={count*step:.17g}
UNST_CFL_NUMBER=0
MGLEVEL=0
FLUID_MODEL=IDEAL_GAS
GAMMA_VALUE=1.4
REF_DIMENSIONALIZATION=DIMENSIONAL
INIT_OPTION=TD_CONDITIONS
FREESTREAM_OPTION=DENSITY_FS
FREESTREAM_DENSITY=1
FREESTREAM_PRESSURE=1
MACH_NUMBER=0.4
PRIMITIVE_REPRESENTATION=NODE
PRIMITIVE_FACE_QUADRATURE=AGGREGATED
PRIMITIVE_VOLUME_QUADRATURE=NODE_POINT
PRIMITIVE_INITIALIZATION={mode}
PRIMITIVE_BOUNDARY_MODE={mode}
PRIMITIVE_INITIAL_STATE=(2,0.5,0.3,3)
PRIMITIVE_BOUNDARY_STATE=(2,0.5,0.3,3)
MARKER_CUSTOM=(boundary)
CONV_NUM_METHOD_FLOW=ROE
MUSCL_FLOW=YES
MUSCL_KAPPA_FLOW=0
NUM_METHOD_GRAD_RECON=WEIGHTED_LEAST_SQUARES
SLOPE_LIMITER_FLOW=NONE
OUTPUT_FILES=(PARAVIEW_MULTIBLOCK,RESTART)
VOLUME_OUTPUT=(COORDINATES,SOLUTION,PRIMITIVE)
VOLUME_FILENAME=primitive
RESTART_FILENAME=primitive_state
CONV_FILENAME=history
OUTPUT_WRT_FREQ={max(1,count//2)}
SCREEN_WRT_FREQ_TIME={max(1,count//2)}
HISTORY_WRT_FREQ_TIME=1
HISTORY_OUTPUT=(ITER,TIME_DOMAIN,PRIMITIVE_RHS,PRIMITIVE_STATE)
OUTPUT_PRECISION=17
PRIMITIVE_GEOMETRY_AUDIT=geometry.csv
'''
    if prefix:
        text=text.replace('VOLUME_FILENAME=primitive','VOLUME_FILENAME='+prefix).replace(
            'RESTART_FILENAME=primitive_state','RESTART_FILENAME='+prefix+'_state').replace(
            'CONV_FILENAME=history','CONV_FILENAME=history_'+prefix).replace(
            'PRIMITIVE_GEOMETRY_AUDIT=geometry.csv','PRIMITIVE_GEOMETRY_AUDIT='+prefix+'_geometry.csv')
    return text


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--n',type=int,default=16,choices=[8,16,32,64])
    a=p.parse_args()
    for directory,family,uniform,stem in [('uniform','regular',True,'uniform'),
            ('vortex','regular',False,'vortex_regular'),('vortex','jitter',False,'vortex_jitter')]:
        path=a.output_dir/directory
        path.mkdir(parents=True,exist_ok=True)
        (path/(stem+'.su2')).write_text(mesh(family,a.n))
        (path/(stem+'.cfg')).write_text(config(stem+'.su2',a.n,uniform,prefix=stem))
    print(a.output_dir.resolve())


if __name__=='__main__':main()
