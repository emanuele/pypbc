"""Basic functions to read and write TrackVis .trk files.
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
    """ Read and parse file header structure.
    """
    header = {}
    for field_name, count, dtype in trk_header_structure:
        header[field_name] = N.fromfile(f, dtype=dtype, count=count)
        #if header[field_name].size==1:
        #    header[field_name] = header[field_name][0]
        #    pass
        # print field_name, ':', header[field_name]
        pass
    assert(f.tell()==1000) # header is always 1000 bytes.
    return header


def write_header(f, header):
    """Write header to file.
    """
    for field_name, count, dtype in trk_header_structure:
        header[field_name].tofile(f)
        pass
    assert(f.tell()==1000) # header is always 1000 bytes.
    return


def print_header(header):
    """Print relevant info of header.
    """
    print "Header:"
    relevant_fields = ['dim', 'voxel_size', 'origin', 'n_count' ]
    for field in relevant_fields:
        print '\t',field, ':', header[field]
        pass
    return


def read_fibers(f, header):
    """Read fibers from file and fill a dictionary.
    """
    fiber = {}
    # structure of each entry of the dictionary:
    # fiber_id : [[X1,Y1,Z1,SCALAR1...],...,[Xn,Yn,Zn,SCALARn...], [PROPERTIES]
    # Note that in PBC2009 trckvis files there are no scalars or
    # properties, which means that the actual structure of the fiber
    # dictionary is simply:
    # fiber_id : [[X1,Y1,Z1],...,[Xn,Yn,Zn]], []
    n_scalars = header['n_scalars'][0]
    n_fibers = header['n_count'][0]
    for fiber_id in range(n_fibers):
        num_points = N.fromfile(f, dtype='<i4', count=1)[0]
        xyz_scalar = N.fromfile(f, dtype='<f4', count=num_points*(3+n_scalars)).reshape(num_points, 3+n_scalars)
        properties = N.fromfile(f, dtype='<f4', count=header['n_properties'][0])
        fiber[fiber_id] = [xyz_scalar, properties]
        if fiber_id%(int(n_fibers/10))==0:
            print 'Reading fibers...', str(1+int(100.0*fiber_id/n_fibers))+'%'
            sys.stdout.flush()
            pass
        pass
    return fiber


def write_fibers(f, fiber):
    n_scalars = header['n_scalars'][0]
    n_fibers = header['n_count'][0]
    for fiber_id in range(n_fibers):
        num_points = N.array((fiber[fiber_id][0]).shape[0], dtype='<i4')
        num_points.tofile(f)
        xyz_scalar = N.array(fiber[fiber_id][0], dtype='<f4')
        xyz_scalar.tofile(f)
        properties = N.array(fiber[fiber_id][1], dtype='<f4')
        properties.tofile(f)
        if fiber_id%(int(n_fibers/10))==0:
            print 'Writing fibers...', str(1+int(100.0*fiber_id/n_fibers))+'%'
            sys.stdout.flush()
            pass
        pass
    return


def mm2voxel(xyz, header):
    """Converts coordinates from mm to voxel.
    """
    return N.floor(xyz/header['voxel_size'])


def voxel2mm(Vxyz, header):
    """Converts coordinates from voxel to mm.
    """
    return (Vxyz+0.5)*header['voxel_size']
    

if __name__=="__main__":
    
    print "This simple program reads a TrackVis .trk file, parse it, build"
    print "structures to represent fibers as Python dictionay of arrays"
    print "and then saves structures in TrackVis .trk file format."
    print "The resulting file is expected to be identical to the original."

    # filename = "dsi.trk"
    # filename = "dti.trk"
    filename = "hardiO10.trk"

    print "file:", filename

    f = open(filename)
    header = read_header(f)
    print_header(header)
    fiber = read_fibers(f, header)
    f.close()
    
    print "Fiber ID=1000:"
    print fiber[1000]
    print "Convert points from mm to voxel coordinate:"
    Vxyz = mm2voxel(fiber[1000][0], header)
    print Vxyz
    print "Convert back and check whether differences are less than grid size...",
    assert(((voxel2mm(Vxyz, header)-fiber[1000][0])<header['voxel_size']).all())
    print "OK."
    
    filename = filename+"_copy"
    print "Saving to:", filename
    f = open(filename,'w')
    write_header(f, header)
    write_fibers(f, fiber)
    f.close()

