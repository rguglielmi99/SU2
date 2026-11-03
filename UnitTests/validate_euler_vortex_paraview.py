#!/usr/bin/env pvpython
"""Read actual diagnostic VTU and PVD with ParaView/VTK, without rendering."""
import argparse
import csv
from pathlib import Path
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.util.numpy_support import vtk_to_numpy
from paraview.simple import PVDReader, Delete
from paraview import servermanager
from euler_vortex_output import point_fields


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vtu-dir',required=True,type=Path)
    p.add_argument('--input-dir',required=True,type=Path)
    p.add_argument('--divisor',type=int,default=8)
    p.add_argument('--report',required=True,type=Path)
    p.add_argument('--errors',type=Path,help='Also independently recompute global norms from VTK error arrays')
    a=p.parse_args()
    records=[]
    errors={}
    if a.errors:
        with a.errors.open() as f:
            errors={(r['family'],int(r['N']),float(r['time']),r['kind'],int(r['component'])):r
                    for r in csv.DictReader(f) if r['group']=='all' and r['variant']==f'FE_d{a.divisor}'}
    for family in ('regular','jitter'):
        for n in (8,16,32,64):
            path=a.input_dir/f'{family}_{n}_{a.divisor}.csv'
            cells=np.loadtxt(str(path)+'.triangles',dtype=int)
            for index,t in enumerate((0.,.25,.5)):
                file=a.vtu_dir/f'{family}_{n}_{index}.vtu'
                reader=vtkXMLUnstructuredGridReader()
                if not reader.CanReadFile(str(file)):raise ValueError('VTK cannot read '+str(file))
                reader.SetFileName(str(file));reader.Update()
                grid=reader.GetOutput()
                data=np.loadtxt(str(path)+f'.snapshot_{t:.6f}')
                assert grid.GetNumberOfPoints()==len(data) and grid.GetNumberOfCells()==len(cells)
                np.testing.assert_array_equal(vtk_to_numpy(grid.GetPoints().GetData())[:,:2],data[:,2:4])
                for name,expected in point_fields(data).items():
                    array=grid.GetPointData().GetArray(name)
                    if array is None:raise ValueError('Missing '+name)
                    np.testing.assert_array_equal(vtk_to_numpy(array),expected)
                for cell in range(len(cells)):
                    ids=grid.GetCell(cell).GetPointIds()
                    assert grid.GetCellType(cell)==5
                    np.testing.assert_array_equal([ids.GetId(k) for k in range(3)],cells[cell])
                assert vtk_to_numpy(grid.GetFieldData().GetArray('TimeValue'))[0]==t
                maximum_norm_difference=0
                if a.errors:
                    v=vtk_to_numpy(grid.GetPointData().GetArray('dual_volume'))
                    for kind,prefix in [('W_primitive','error_{}_primitive'),('W_conservative','error_{}_conservative'),
                                        ('difference_W','difference_{}_primitive_conservative')]:
                        for k,name in enumerate(('rho','u','v','p')):
                            e=vtk_to_numpy(grid.GetPointData().GetArray(prefix.format(name)))
                            measured=dict(rms=np.sqrt(np.mean(e*e)),l1=np.sum(v*abs(e))/v.sum(),
                                          l2=np.sqrt(np.sum(v*e*e)/v.sum()),linf=np.max(abs(e)))
                            reference=errors[family,n,t,kind,k]
                            for norm,value in measured.items():
                                expected=float(reference[norm])
                                np.testing.assert_allclose(value,expected,rtol=2e-14,atol=2e-16)
                                maximum_norm_difference=max(maximum_norm_difference,abs(value-expected))
                records.append(dict(family=family,N=n,time=t,reader='VTK XML',points=grid.GetNumberOfPoints(),
                                    cells=grid.GetNumberOfCells(),values='exact',connectivity='exact',
                                    norms_verified=int(bool(a.errors)),max_norm_difference=maximum_norm_difference))
            sequence=PVDReader(FileName=str(a.vtu_dir/f'{family}_{n}.pvd'))
            assert list(sequence.TimestepValues)==[0.,.25,.5]
            for t in (0.,.25,.5):
                sequence.UpdatePipeline(time=t)
                output=servermanager.Fetch(sequence)
                # PVD exposes a multiblock in some ParaView versions.
                if output.IsA('vtkMultiBlockDataSet'):output=output.GetBlock(0)
                exact=np.loadtxt(str(path)+f'.snapshot_{t:.6f}')
                np.testing.assert_array_equal(vtk_to_numpy(output.GetPointData().GetArray('p_exact')),exact[:,23])
            Delete(sequence)
    with a.report.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    print('Actual ParaView/VTK reading passed: 24 VTU files and 8 time collections; no GUI/rendering claimed.')


if __name__=='__main__':main()
