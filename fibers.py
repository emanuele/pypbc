"""Basic functions to read and write TrackVis .trk files and to play
with fibers.

Copyright (c) 2009 Emanuele Olivetti <emanuele_AT_relativita.com>

This library is free software; you can redistribute it and/or modify
it either under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.
"""

import numpy as N
import sys

# Definition of trackvis header structure.
# See http://www.trackvis.org/docs/?subsect=fileformat
# See http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html
trk_header_structure = [['id_string', 1, 'S6'],
                        ['dim', 3, '<h'],
                        ['voxel_size', 3, '<f4'],
                        ['origin', 3, '<f4'],
                        ['n_scalars', 1, '<h'],
                        ['scalar_name', 10, 'S20'],
                        ['n_properties', 1, '<h'],
                        ['property_name', 10, 'S20'],
                        ['reserved', 1, 'S508'],
                        ['voxel_order', 1, 'S4'],
                        ['pad2', 1, 'S4'],
                        ['image_orientation_patient', 6, '<f4'],
                        ['pad1', 1, 'S2'],
                        ['invert_x', 1, 'S1'],
                        ['invert_y', 1, 'S1'],
                        ['invert_z', 1, 'S1'],
                        ['swap_xy', 1, 'S1'],
                        ['swap_yz', 1, 'S1'],
                        ['swap_zx', 1, 'S1'],
                        ['n_count', 1, '<i4'],
                        ['version', 1, '<i4'],
                        ['hdr_size', 1, '<i4'],
                        ]


def read_header(f):
    """ Read and parse .trk file header structure.
    """
    header = {}
    for field_name, count, dtype in trk_header_structure:
        header[field_name] = N.fromfile(f, dtype=dtype, count=count)
        pass
    assert(f.tell()==1000) # header is always 1000 bytes.
    return header


def write_header(f, header):
    """Write .trk header to file.
    """
    for field_name, count, dtype in trk_header_structure:
        # Note that ".astype(dtype)" is just to be sure or correct types:
        header[field_name].astype(dtype).tofile(f)
        pass
    assert(f.tell()==1000) # header is always 1000 bytes.
    return


def print_header(header):
    """Print relevant info of .trk header.
    """
    print "Header:"
    relevant_fields = ['dim', 'voxel_size', 'origin', 'n_count' ]
    for field in relevant_fields:
        print '\t',field, ':', header[field]
        pass
    return


def progress_meter(position, total, message, steps=10):
    """Simple progress meter.
    """
    if position%(int(total/steps))==0:
        print message, str(1+int(100.0*position/total))+'%'
        sys.stdout.flush()
        pass
    return


def read_fibers(f, header):
    """Read fibers from .trk file and fill a list.
    """
    fiber = []
    # structure of each entry of the list:
    # [[X1,Y1,Z1,SCALAR1...],...,[Xn,Yn,Zn,SCALARn...]], [PROPERTIES]
    # Note that in PBC2009 trckvis files there are no scalars or
    # properties, which means that the actual structure of the fiber
    # list is simply:
    # fiber_id : [[X1,Y1,Z1],...,[Xn,Yn,Zn]], []
    n_scalars = header['n_scalars'][0]
    n_fibers = header['n_count'][0]
    for fiber_id in range(n_fibers):
        num_points = N.fromfile(f, dtype='<i4', count=1)[0]
        xyz_scalar = N.fromfile(f, dtype='<f4', count=num_points*(3+n_scalars)).reshape(num_points, 3+n_scalars)
        properties = N.fromfile(f, dtype='<f4', count=header['n_properties'][0])
        fiber.append([xyz_scalar, properties])
        progress_meter(fiber_id, n_fibers, 'Reading fibers...')
        pass
    return fiber


def write_fibers(f, fiber, header):
    """Write fibers to file in .trk format. Assumption: header has
    already been written.
    """
    n_scalars = header['n_scalars'][0]
    n_fibers = header['n_count'][0]
    for fiber_id in range(n_fibers):
        num_points = N.array((fiber[fiber_id][0]).shape[0], dtype='<i4')
        num_points.tofile(f)
        xyz_scalar = N.array(fiber[fiber_id][0], dtype='<f4')
        xyz_scalar.tofile(f)
        properties = N.array(fiber[fiber_id][1], dtype='<f4')
        properties.tofile(f)
        progress_meter(fiber_id, n_fibers, 'Writing fibers...')
        pass
    return


def mm2voxel(xyz, header):
    """Converts coordinates from mm to voxel.
    """
    return N.floor(xyz/header['voxel_size']).astype('i')


def voxel2mm(Vxyz, header):
    """Converts coordinates from voxel to mm.
    """
    return (Vxyz+0.5)*header['voxel_size']
    

def build_voxel_fibers_dict(fiber, header):
    """Build a dictionary that given a voxel returns all fibers (IDs)
    crossing it.
    """
    voxel2fibers = {}
    n_fibers = len(fiber)
    for fiber_id in range(n_fibers):
        xyz = fiber[fiber_id][0]
        ijk = mm2voxel(xyz, header)
        for i in range(xyz.shape[0]):
            try:
                voxel2fibers[tuple(ijk[i,:])].append(fiber_id)
            except KeyError:
                voxel2fibers[tuple(ijk[i,:])] = [fiber_id]
                pass
            pass
        progress_meter(fiber_id, n_fibers, 'Mapping voxels to fibers...')
        pass
    n_voxels = len(voxel2fibers.keys())
    # Now transform each list of IDs in an array of IDs:
    for n, ijk in enumerate(voxel2fibers.keys()):
        voxel2fibers[ijk] = N.array(voxel2fibers[ijk])
        progress_meter(n, n_voxels, 'Converting lists to arrays...')
        pass
    return voxel2fibers


if __name__=="__main__":
    
    print "This simple program reads a TrackVis .trk file, parse it, build"
    print "structures to represent fibers as Python list of arrays"
    print "and then saves structures in TrackVis .trk file format."
    print "The resulting file is expected to be identical to the original."
    print "As a further step a dictionary, mapping voxel to fibers, is built"
    print "and some examples using it are shown."
    
    # filename = "dsi.trk"
    # filename = "dti.trk"
    filename = "hardiO10.trk"
    
    print
    print "file:", filename
    f = open(filename)
    header = read_header(f)
    print_header(header)
    fiber = read_fibers(f, header)
    f.close()
    
    print
    fiber_id = 1000
    print "Example: fiber_id=",fiber_id
    print fiber[fiber_id]
    print "Convert points from mm to voxel coordinates:"
    Vxyz = mm2voxel(fiber[fiber_id][0], header)
    print Vxyz
    print "Convert back and check whether differences are less than grid size...",
    assert(((voxel2mm(Vxyz, header)-fiber[fiber_id][0])<header['voxel_size']).all())
    print "OK."
    
    print
    filename2 = filename+"_COPY.trk"
    print "Saving to:", filename2
    f = open(filename2,'w')
    write_header(f, header)
    write_fibers(f, fiber, header)
    f.close()
    
    print
    print "Building voxel2fibers dictionary:"
    voxel2fibers = build_voxel_fibers_dict(fiber, header)
    voxel = tuple(header['dim'] / 2)
    print "Example: fibers crossing voxel", voxel
    print voxel2fibers[voxel]
    
    print
    x = header['dim'][0] / 2
    print "Example: counting fibers crossing plane x =", x
    counter = 0
    for y in range(header['dim'][1]):
        for z in range(header['dim'][2]):
            try:
                counter += voxel2fibers[(x,y,z)].size
            except KeyError:
                pass
            pass
        pass
    print "Number of fibers:", counter
    
    print
    fiber_id = 2000
    print "Which fibers cross (the voxels of) fiber[fiber_id=",fiber_id,"] ?"
    xyz = fiber[fiber_id][0]
    ijk = mm2voxel(xyz, header)
    fiber_id_list = N.unique(N.hstack([voxel2fibers[i,j,k] for i,j,k in ijk]))
    print fiber_id_list
    print fiber_id_list.size, "fibers."

    print
    print "Saving .trk file with just the previous list of fibers."
    filename3 = filename+'_cross_fiber_id_'+str(fiber_id)+'.trk'
    print "Saving to:", filename3
    import copy
    fiber2 = [fiber[fiber_id] for fiber_id in fiber_id_list]
    header2 = copy.deepcopy(header)
    header2['n_count'] = N.array([fiber_id_list.size])
    f = open(filename3, 'w')
    write_header(f, header2)
    write_fibers(f, fiber2, header2)
    f.close()
