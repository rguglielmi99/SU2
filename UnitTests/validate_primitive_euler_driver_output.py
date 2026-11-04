#!/usr/bin/env pvpython
"""Validate actual SU2 Float32 VTU against native Float64 checkpoints/history."""
import argparse
import csv
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader, vtkXMLMultiBlockDataReader
from vtkmodules.util.numpy_support import vtk_to_numpy
from run_primitive_euler_driver_e2e import binary, reference


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',required=True,type=Path)
    p.add_argument('--report',required=True,type=Path)
    a=p.parse_args();results=[]
    for case in sorted(a.input_dir.glob('*_*_*')):
        if not (case/'case.cfg').exists() or case.name.startswith('negative_'):continue
        kind,family,n=case.name.split('_');n=int(n)
        if kind not in ('vortex','uniform'):continue
        with (case/'history.csv').open() as f:
            reader=csv.DictReader(f,skipinitialspace=True);reader.fieldnames=[k.strip() for k in reader.fieldnames]
            history={int(float(r['Time_Iter'])):r for r in reader}
        with (case/'geometry.csv').open() as f:
            volume={int(r['first']):float(r['volume']) for r in csv.DictReader(f) if r['kind']=='node'}
        for vtm in sorted(case.glob('primitive_*.vtm')):
            index=int(vtm.stem.rsplit('_',1)[1]);t=float(history[index]['Cur_Time'])
            xml=ET.parse(vtm).getroot()
            meta=float(xml.find('FieldData/DataArray').text)
            np.testing.assert_allclose(meta,t,rtol=0,atol=1e-14)
            dataset=next(d for d in xml.findall('.//DataSet') if d.attrib['name']=='Internal')
            vtu=case/dataset.attrib['file']
            reader=vtkXMLUnstructuredGridReader();reader.SetFileName(str(vtu));reader.Update();g=reader.GetOutput()
            assert g.GetNumberOfPoints()==(n+1)**2 and g.GetNumberOfCells()==2*n*n
            b=binary(case/f'primitive_state_{index:05d}.dat')
            ids=vtk_to_numpy(g.GetPointData().GetArray('Global_ID')).astype(int)
            np.testing.assert_array_equal(ids,np.arange((n+1)**2))
            points=vtk_to_numpy(g.GetPoints().GetData())[:,:2]
            exactpoints=np.column_stack((b['x'],b['y']))
            np.testing.assert_array_equal(points,exactpoints.astype(np.float32))
            fields={name:vtk_to_numpy(g.GetPointData().GetArray(name)) for name in
                    ('Density','Velocity','Pressure','Dual_Volume','Physical_Time','Entropy','Temperature','Speed','Mach')}
            state=np.column_stack((fields['Density'],fields['Velocity'][:,0],fields['Velocity'][:,1],fields['Pressure']))
            exact=np.column_stack([b[k] for k in ('Density','Velocity_x','Velocity_y','Pressure')])
            np.testing.assert_array_equal(state,exact.astype(np.float32))
            np.testing.assert_array_equal(fields['Physical_Time'],np.full(len(state),np.float32(t)))
            np.testing.assert_array_equal(fields['Dual_Volume'],np.array([volume[i] for i in ids],dtype=np.float32))
            assert all(np.all(np.isfinite(f)) for f in fields.values())
            assert np.min(fields['Temperature'])>0
            for cell in range(g.GetNumberOfCells()):
                assert g.GetCellType(cell)==5
            # Compare triangle sets in global-ID space, allowing only ordering changes.
            tri=[]
            for cell in range(g.GetNumberOfCells()):
                local=g.GetCell(cell).GetPointIds();tri.append(tuple(sorted(ids[local.GetId(k)] for k in range(3))))
            archived=ET.parse(Path(__file__).resolve().parents[1]/f'Docs/primitive_euler_phase7_vortex/{family}_{n}_0.vtu').getroot()
            c=next(x for x in archived.findall('.//Cells/DataArray') if x.attrib['Name']=='connectivity')
            expected=np.fromstring(c.text,sep=' ',dtype=int).reshape(-1,3)
            assert sorted(tri)==sorted(tuple(sorted(row)) for row in expected)
            # Metadata and block references are also read by the actual VTK multiblock reader.
            multi=vtkXMLMultiBlockDataReader();multi.SetFileName(str(vtm));multi.Update()
            assert multi.GetOutput().GetNumberOfBlocks()>0
            quantization=np.max(np.abs(state.astype(float)-exact),axis=0)
            if kind=='vortex':
                archive_index={0.:0,.25:1,.5:2}[t]
                _,ref=reference(family,n,archive_index)
                analytic=np.column_stack([ref[k+'_exact'] for k in ('rho','u','v','p')])
            else:analytic=np.tile((2,.5,.3,3),(len(state),1))
            weights=np.array([volume[i] for i in ids])
            l2_double=np.sqrt(np.sum(weights[:,None]*(exact-analytic)**2,axis=0)/weights.sum())
            l2_vtk=np.sqrt(np.sum(weights[:,None]*(state.astype(float)-analytic)**2,axis=0)/weights.sum())
            bound=np.sqrt(np.sum(weights[:,None]*(state.astype(float)-exact)**2,axis=0)/weights.sum())
            if np.max(abs(l2_double-l2_vtk)-bound)>2e-15:
                raise ValueError('VTU norms differ beyond their measured quantization bound')
            for k in range(4):
                results.append(dict(case=case.name,time=t,component=k,max_state_quantization=quantization[k],
                    points=g.GetNumberOfPoints(),triangles=g.GetNumberOfCells(),vtk_values_equal_float32_cast=1,
                    geometry_equal_float32_cast=1,time_equal_float32_cast=1,l2_from_checkpoint=l2_double[k],
                    l2_from_vtu=l2_vtk[k],l2_quantization_bound=bound[k],l2_norm_difference=abs(l2_double[k]-l2_vtk[k])))
    if not results:raise ValueError('No native outputs read')
    with a.report.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    print('Native VTM/VTU read and validated:',len(results)//4,'snapshots; exact Float32 casts of Float64 state.')


if __name__=='__main__':main()
