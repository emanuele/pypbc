# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""Basic functions to read and write TrackVis .trk files and to play
with streamlines.

Copyright (c) 2009 Emanuele Olivetti <emanuele_AT_relativita.com>

This library is free software; you can redistribute it and/or modify
it either under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.
"""

import numpy as N
import sys
import copy
import operator

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



def progress_meter(position, total, message, steps=10):
    """Simple progress meter.
    """
    if position%(int(total/steps))==0:
        print message, str(1+int(100.0*position/total))+'%'
        sys.stdout.flush()
        pass
    return


class Streamlines(object):
    """Class to deal with some streamlines.
    """
    def __init__(self, header=None, streamline=None):
        """Constructor.
        """
        self.header = header
        if self.header == None:
            self.header = {}
            pass
        self.streamline = streamline
        if self.streamline == None:
            self.streamline = {}
            pass
        self.filename = None
        self.voxel2streamlines = None
        return

    def loadTrk(self, filename):
        """Load streamlines from TrackVis .trk file.
        """
        self.filename = filename
        f = open(filename)
        self.readHeaderTrk(f)
        self.readStreamlinesTrk(f)
        f.close()
        return self.streamline

    def saveTrk(self, filename):
        """Save streamlines to a TrackVis .trk file.
        """
        self.filename = filename
        f = open(filename, 'w')
        self.writeHeaderTrk(f)
        self.writeStreamlinesTrk(f)
        f.close()
        return

    def readHeaderTrk(self, f):
        """Read and parse .trk file header structure.
        """
        global trk_header_structure
        self.header = {}
        for field_name, count, dtype in trk_header_structure:
            self.header[field_name] = N.fromfile(f, dtype=dtype, count=count)
            pass
        assert(f.tell()==1000) # header is always 1000 bytes.
        return

    def printHeaderTrk(self):
        """Print relevant info of .trk header.
        """
        print "Header:"
        relevant_fields = ['dim', 'voxel_size', 'origin', 'n_count' ]
        for field in relevant_fields:
            print '\t',field, ':', self.header[field]
            pass
        return

    def writeHeaderTrk(self, f):
        """Write .trk header to file.
        """
        global trk_header_structure
        for field_name, count, dtype in trk_header_structure:
            # Note that ".astype(dtype)" is just to be sure or correct types:
            self.header[field_name].astype(dtype).tofile(f)
            pass
        assert(f.tell()==1000) # header is always 1000 bytes.
        return

    def readStreamlinesTrk(self, f):
        """Read streamlines from .trk file and fill a list.
        """
        self.streamline = []
        # structure of each entry of the self.streamline list:
        # [[X1,Y1,Z1,SCALAR1...],...,[Xn,Yn,Zn,SCALARn...]], [PROPERTIES]
        # Note that in PBC2009 trckvis files there are no scalars or
        # properties, which means that the actual structure of the streamline
        # list is simply:
        # streamline_id : [[X1,Y1,Z1],...,[Xn,Yn,Zn]], []
        n_scalars = self.header['n_scalars'][0]
        n_streamlines = self.header['n_count'][0]
        for streamline_id in range(n_streamlines):
            num_points = N.fromfile(f, dtype='<i4', count=1)[0]
            xyz_scalar = N.fromfile(f, dtype='<f4', count=num_points*(3+n_scalars)).reshape(num_points, 3+n_scalars)
            properties = N.fromfile(f, dtype='<f4', count=self.header['n_properties'][0])
            self.streamline.append([xyz_scalar, properties])
            progress_meter(streamline_id, n_streamlines, 'Reading streamlines...')
            pass
        return
        
    def writeStreamlinesTrk(self, f):
        """Write streamlines to file in .trk format. Assumption: header has
        already been written.
        """
        n_scalars = self.header['n_scalars'][0]
        n_streamlines = self.header['n_count'][0]
        for streamline_id in range(n_streamlines):
            num_points = N.array((self.streamline[streamline_id][0]).shape[0], dtype='<i4')
            num_points.tofile(f)
            xyz_scalar = N.array(self.streamline[streamline_id][0], dtype='<f4')
            xyz_scalar.tofile(f)
            properties = N.array(self.streamline[streamline_id][1], dtype='<f4')
            properties.tofile(f)
            progress_meter(streamline_id, n_streamlines, 'Writing streamlines...')
            pass
        return

    def mm2voxel(self, xyz):
        """Converts coordinates from mm to voxel.
        """
        return N.floor(xyz/self.header['voxel_size']).astype('i')

    def voxel2mm(self, Vxyz):
        """Converts coordinates from voxel to mm.
        """
        return (Vxyz+0.5)*self.header['voxel_size']

    def buildVoxelStreamlinesDict(self):
        """Build a dictionary that given a voxel returns all streamlines (IDs)
        crossing it.
        """
        self.voxel2streamlines = {}
        n_streamlines = len(self.streamline)
        for streamline_id in range(n_streamlines):
            xyz = self.streamline[streamline_id][0]
            ijk = self.mm2voxel(xyz)
            for i in range(xyz.shape[0]):
                try:
                    self.voxel2streamlines[tuple(ijk[i,:])].append(streamline_id)
                except KeyError:
                    self.voxel2streamlines[tuple(ijk[i,:])] = [streamline_id]
                    pass
                pass
            progress_meter(streamline_id, n_streamlines, 'Mapping voxels to streamlines...')
            pass
        n_voxels = len(self.voxel2streamlines.keys())
        # Now transform each list of IDs in an array of IDs:
        for n, ijk in enumerate(self.voxel2streamlines.keys()):
            self.voxel2streamlines[ijk] = N.array(self.voxel2streamlines[ijk])
            progress_meter(n, n_voxels, 'Converting lists to arrays...')
            pass
        return self.voxel2streamlines

    def fastCopy(self):
        """Create a copy of the current object avoiding the self.fiber
        list, which could be quite expensive.
        """
        new_streamlines = Streamlines()
        new_streamlines.header = copy.deepcopy(self.header)
        return new_streamlines

    def selectStreamlines(self, streamlines_ids):
        """Given a list of streamlines ID returns a new Streamline
        object with just those streamlines and an header equivalent to
        the current one except for header['n_count'], which is
        adjusted to the actual number of fibers.
        """
        new_streamlines = self.fastCopy()
        new_streamlines.streamline = [self.streamline[streamline_id] for streamline_id in streamlines_ids]
        new_streamlines.header['n_count'] = N.array([len(new_streamlines.streamline)]).astype('<i4')
        return new_streamlines

    def selectStreamlinesFromVoxels(self, voxels):
        """Select streamlines specifying a list of voxels they must cross.
        """
        # WARNING: if a voxels has no fibers than KeyError is raised!
        #
        # SLOW:
        # tmp = [list(self.voxel2streamlines[tuple(voxel)]) for voxel in voxels]
        # streamlines_ids = unique(reduce(operator.add, tmp))
        #
        # FAST:
        streamlines_ids = N.unique(N.hstack([self.voxel2streamlines[tuple(v)] for v in voxels]))
        # ALTERNATIVE:
        # streamlines_ids = N.unique(N.hstack([self.voxel2streamlines[i,j,k] for i,j,k in voxels]))
        return self.selectStreamlines(streamlines_ids)

    def getVolume(self):
        """Return a volume where a voxel is 1 if at least one fiber
        crosses it, otherwise 0.
        """
        volume = N.zeros(self.header['dim'], dtype='i')
        if self.voxel2streamlines is not None: # fast:
            ijk = N.array(self.voxel2streamlines.keys())
            volume[ijk[:,0],ijk[:,1],ijk[:,2]] = 1.0
            pass
        else: # slow but does not require voxel2streamlines:
            for xyz, tmp in self.streamline:
                ijk = self.mm2voxel(xyz)
                volume[ijk[:,0],ijk[:,1],ijk[:,2]] = 1
                pass
            pass
        return volume


if __name__=="__main__":
    
    pass
