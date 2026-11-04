#!/usr/bin/env python3
"""Generate native physical-BC examples on the archived regular/jitter mesh."""
import argparse
from pathlib import Path
import math
from primitive_euler_driver_cases import mesh, config


def case(name, family='regular', n=16):
    text=mesh(family,n)
    head, boundary=text.split('NMARK= 1\n')
    edges=[line for line in boundary.splitlines() if line.startswith('3 ')]
    wall=[]; inlet=[]; outlet=[]
    for line in edges:
        i,j=map(int,line.split()[1:])
        if i//(n+1)==j//(n+1) and i//(n+1) in (0,n): wall.append(line)
        elif i%(n+1)==j%(n+1)==0: inlet.append(line)
        else: outlet.append(line)
    supersonic='supersonic' in name
    speed=2.0 if supersonic else .4
    mode='sym' if name=='symmetry' else ('wall' if name in ('wall','walls_farfield','walls_inlet_outlet') else 'open')
    groups={'lateral':wall,'inlet':inlet,'outlet':outlet}
    mesh_text=head+f'NMARK= {len(groups)}\n'
    for tag,lines in groups.items():
        mesh_text+=f'MARKER_TAG= {tag}\nMARKER_ELEMS= {len(lines)}\n'+'\n'.join(lines)+'\n'
    cfg=config(name+'.su2',n,True,steps=20,dt=.001,prefix=name)
    cfg=cfg.replace('% Native single-zone primitive solver; no source, no general conservative BC.',
                    '% Physical Euler BCs: native marker names and dimensional ideal-gas data.')
    cfg='\n'.join(line for line in cfg.splitlines() if not line.startswith(('MARKER_CUSTOM=', 'PRIMITIVE_BOUNDARY_MODE=', 'PRIMITIVE_BOUNDARY_STATE=')))+'\n'
    cfg=cfg.replace('PRIMITIVE_INITIAL_STATE=(2,0.5,0.3,3)\n','').replace('PRIMITIVE_INITIALIZATION=UNIFORM','PRIMITIVE_INITIALIZATION=FREESTREAM')
    cfg=cfg.replace('PRIMITIVE_INITIALIZATION=FREESTREAM\n',
                    'PRIMITIVE_INITIALIZATION=FREESTREAM\nPRIMITIVE_BOUNDARY_MODE=LOCAL_CHARACTERISTIC\n')
    cfg=cfg.replace('MACH_NUMBER=0.4',f'MACH_NUMBER={speed/math.sqrt(1.4):.17g}')
    cfg+='GAS_CONSTANT=1\n'
    if name.startswith('farfield'): cfg+='MARKER_FAR=(lateral,inlet,outlet)\n'
    else:
        cfg+=('MARKER_SYM=(lateral)\n' if mode=='sym' else ('MARKER_EULER=(lateral)\n' if mode=='wall' else 'MARKER_FAR=(lateral)\n'))
        if name in ('wall','symmetry','walls_farfield'): cfg+='MARKER_FAR=(inlet,outlet)\n'
        elif supersonic:
            cfg+=f'MARKER_SUPERSONIC_INLET=(inlet,1,1,{speed},0,0)\nMARKER_SUPERSONIC_OUTLET=(outlet)\n'
        else:
            tt=1+.4*speed*speed/(2*1.4);pt=tt**3.5
            cfg+=f'INLET_TYPE=TOTAL_CONDITIONS\nMARKER_INLET=(inlet,{tt:.17g},{pt:.17g},1,0,0)\nMARKER_OUTLET=(outlet,1)\n'
    return cfg,mesh_text


NAMES=('wall','symmetry','farfield_subsonic','farfield_supersonic',
       'inlet_outlet_subsonic','inlet_outlet_supersonic','walls_farfield','walls_inlet_outlet')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--family',choices=('regular','jitter'),default='regular')
    parser.add_argument('--n',type=int,choices=(8,16,32,64),default=16)
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    for name in NAMES:
        cfg,grid=case(name,args.family,args.n)
        shared='channel_'+args.family+'.su2'
        (args.output_dir/(name+'.cfg')).write_text(cfg.replace('MESH_FILENAME='+name+'.su2','MESH_FILENAME='+shared))
        (args.output_dir/shared).write_text(grid)


if __name__=='__main__':main()
