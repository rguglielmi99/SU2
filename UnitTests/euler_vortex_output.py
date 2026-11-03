#!/usr/bin/env python3
"""Diagnostic ASCII VTU/PVD writer and strict round-trip validation."""
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np


def point_fields(a):
    result = {'node_id': a[:, 0].astype(int), 'dual_volume': a[:, 1]}
    for label, start in [('primitive', 4), ('conservative', 12), ('exact', 20)]:
        for k, name in enumerate(('rho', 'u', 'v', 'p')):
            result[f'{name}_{label}'] = a[:, start+k]
        result[f'velocity_magnitude_{label}'] = np.hypot(a[:, start+1], a[:, start+2])
        result[f'entropy_{label}'] = a[:, start+3]/a[:, start]**1.4
    for k, name in enumerate(('rho', 'u', 'v', 'p')):
        result[f'error_{name}_primitive'] = a[:, 4+k]-a[:, 20+k]
        result[f'error_{name}_conservative'] = a[:, 12+k]-a[:, 20+k]
        result[f'difference_{name}_primitive_conservative'] = a[:, 4+k]-a[:, 12+k]
    if not all(np.all(np.isfinite(x)) for x in result.values()):
        raise ValueError('Nonfinite export')
    return result


def array(parent, name, values, components=1, integer=False):
    a = ET.SubElement(parent, 'DataArray', type='Int64' if integer else 'Float64', Name=name,
                      NumberOfComponents=str(components), format='ascii')
    a.text = ' '.join(str(int(x)) if integer else format(float(x), '.17g') for x in np.asarray(values).ravel())
    return a


def write_vtu(path, a, triangles, time):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    root = ET.Element('VTKFile', type='UnstructuredGrid', version='0.1', byte_order='LittleEndian')
    grid = ET.SubElement(root, 'UnstructuredGrid')
    array(ET.SubElement(grid, 'FieldData'), 'TimeValue', [time]).set('NumberOfTuples', '1')
    piece = ET.SubElement(grid, 'Piece', NumberOfPoints=str(len(a)), NumberOfCells=str(len(triangles)))
    pd = ET.SubElement(piece, 'PointData', Scalars='p_primitive')
    for name, values in point_fields(a).items():
        array(pd, name, values, integer=name == 'node_id')
    array(ET.SubElement(piece, 'Points'), 'Points', np.column_stack((a[:, 2:4], np.zeros(len(a)))), 3)
    cells = ET.SubElement(piece, 'Cells')
    array(cells, 'connectivity', triangles, integer=True)
    array(cells, 'offsets', 3*np.arange(1, len(triangles)+1), integer=True)
    array(cells, 'types', np.full(len(triangles), 5), integer=True)
    ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)
    validate_vtu(path, a, triangles, time)


def validate_vtu(path, a, triangles, time):
    root = ET.parse(path).getroot()
    if root.attrib['type'] != 'UnstructuredGrid':
        raise ValueError('Wrong VTK type')
    grid = root.find('UnstructuredGrid')
    piece = grid.find('Piece')
    if int(piece.attrib['NumberOfPoints']) != len(a) or int(piece.attrib['NumberOfCells']) != len(triangles):
        raise ValueError('Wrong VTK dimensions')
    points = np.fromstring(piece.find('Points/DataArray').text, sep=' ').reshape(-1, 3)
    np.testing.assert_array_equal(points[:, :2], a[:, 2:4])
    np.testing.assert_array_equal(points[:, 2], 0)
    cells = {x.attrib['Name']: np.fromstring(x.text, sep=' ', dtype=int) for x in piece.findall('Cells/DataArray')}
    np.testing.assert_array_equal(cells['connectivity'].reshape(-1, 3), triangles)
    np.testing.assert_array_equal(cells['offsets'], 3*np.arange(1, len(triangles)+1))
    np.testing.assert_array_equal(cells['types'], 5)
    if np.min(triangles) < 0 or np.max(triangles) >= len(a):
        raise ValueError('Invalid connectivity')
    exported = {x.attrib['Name']: np.fromstring(x.text, sep=' ') for x in piece.findall('PointData/DataArray')}
    if exported.keys() != point_fields(a).keys():
        raise ValueError('Missing/extra point fields')
    for name, expected in point_fields(a).items():
        np.testing.assert_array_equal(exported[name], expected)
    if float(grid.find('FieldData/DataArray').text) != time:
        raise ValueError('Wrong snapshot time')
    return exported


def write_pvd(path, snapshots):
    root = ET.Element('VTKFile', type='Collection', version='0.1', byte_order='LittleEndian')
    collection = ET.SubElement(root, 'Collection')
    times = []
    for t, file in snapshots:
        if not (Path(path).parent/file).exists():
            raise ValueError('Missing indexed snapshot')
        ET.SubElement(collection, 'DataSet', timestep=str(t), group='', part='0', file=file)
        times.append(t)
    if times != sorted(set(times)):
        raise ValueError('Nonunique or unordered snapshot indices')
    ET.ElementTree(root).write(path, encoding='utf-8', xml_declaration=True)
    parsed = ET.parse(path).getroot().findall('Collection/DataSet')
    assert [(float(x.attrib['timestep']), x.attrib['file']) for x in parsed] == list(snapshots)
