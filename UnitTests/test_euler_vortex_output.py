#!/usr/bin/env python3
"""Writer failures and signed-error semantics, independent of convergence order."""
import tempfile
import unittest
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from euler_vortex_output import write_vtu, validate_vtu, write_pvd, point_fields


class VortexOutput(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.path=Path(self.directory.name)
        self.a=np.zeros((3,28))
        self.a[:,0]=np.arange(3)
        self.a[:,1]=1/6
        self.a[:,2:4]=[[0,0],[1,0],[0,1]]
        self.a[:,4:8]=[[1,1,1,.9],[1.1,1,1,1],[.9,1,1,1.1]]
        self.a[:,12:16]=self.a[:,4:8]+.01
        self.a[:,20:24]=1
        self.cells=np.array([[0,1,2]])

    def tearDown(self):
        self.directory.cleanup()

    def test_signed_fields_and_roundtrip(self):
        write_vtu(self.path/'a.vtu',self.a,self.cells,.25)
        fields=validate_vtu(self.path/'a.vtu',self.a,self.cells,.25)
        self.assertLess(fields['error_p_primitive'][0],0)
        self.assertGreater(fields['error_p_primitive'][2],0)
        np.testing.assert_allclose(fields['difference_p_primitive_conservative'],-.01)
        write_pvd(self.path/'a.pvd',[(.25,'a.vtu')])

    def test_wrong_connectivity_rejected(self):
        write_vtu(self.path/'a.vtu',self.a,self.cells,0)
        tree=ET.parse(self.path/'a.vtu')
        tree.find('.//Cells/DataArray[@Name="connectivity"]').text='0 2 1'
        tree.write(self.path/'a.vtu')
        with self.assertRaises(AssertionError):validate_vtu(self.path/'a.vtu',self.a,self.cells,0)

    def test_wrong_time_and_missing_snapshot_rejected(self):
        write_vtu(self.path/'a.vtu',self.a,self.cells,0)
        with self.assertRaises(ValueError):validate_vtu(self.path/'a.vtu',self.a,self.cells,.5)
        with self.assertRaises(ValueError):write_pvd(self.path/'bad.pvd',[(0,'missing.vtu')])
        with self.assertRaises(ValueError):write_pvd(self.path/'bad.pvd',[(.5,'a.vtu'),(0,'a.vtu')])

    def test_nonfinite_field_rejected(self):
        self.a[0,4]=np.nan
        with self.assertRaises(ValueError):point_fields(self.a)


if __name__=='__main__':unittest.main()
