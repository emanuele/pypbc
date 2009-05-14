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
import gzip
import cPickle
import os

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


class Streamlines(object):
    """Class to deal with some streamlines.
    """

    unknown_streamline = 0
    
    def __init__(self, header=None, streamline=None, properties=None, streamline_id=None, streamlineid_pos=None, streamline_label=None, num_label=None, label_num=None):
        """Constructor.
        """
        self.header = header
        if self.header == None: self.header = {}
        self.streamline = streamline
        if self.streamline == None: self.streamline = []
        self.properties = properties
        if self.properties == None: self.properties = []
        self.streamline_id = streamline_id
        if self.streamline_id == None: self.streamline_id = []
        self.streamlineid_pos = streamlineid_pos
        if self.streamlineid_pos == None: self.streamlineid_pos = []
        self.streamline_label = streamline_label
        if self.streamline_label == None: self.streamline_label = []
        self.num_label = num_label
        if self.num_label == None: self.num_label = {}
        self.label_num = label_num
        if self.label_num == None: self.label_num = {}
        self.filename = None
        self.voxel2streamlines = None
        self.streamline2streamlines = None
        return

    @staticmethod
    def progress_meter(position, total, message, steps=10):
        """Simple progress meter. Static method.
        """
        try:
            if position%(int(total/steps))==0:
                print message, str(1+int(100.0*position/total))+'%'
                sys.stdout.flush()
                pass
        except ZeroDivisionError: # it happens when (int(total/steps))==0
            pass
        return

    def loadTrk(self, filename):
        """Load streamlines from TrackVis .trk file.
        """
        print "Loading",filename,"..."
        self.filename = filename
        f = open(filename)
        self.readHeaderTrk(f)
        self.readStreamlinesTrk(f)
        f.close()
        print "Done."
        return self.streamline

    def saveTrk(self, filename):
        """Save streamlines to a TrackVis .trk file.
        """
        print "Saving",filename,"..."
        self.filename = filename
        f = open(filename, 'w')
        self.writeHeaderTrk(f)
        self.writeStreamlinesTrk(f)
        f.close()
        print "Done."
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
        self.properties = []
        # structure of each entry of the self.streamline list:
        # [[X1,Y1,Z1,SCALAR1...],...,[Xn,Yn,Zn,SCALARn...]]
        n_scalars = self.header['n_scalars'][0]
        n_streamlines = self.header['n_count'][0]
        self.streamline_id = range(1, n_streamlines+1)
        self.streamlineid_pos = dict(zip(self.streamline_id, range(n_streamlines)))
        self.streamline_label = N.ones(n_streamlines, dtype=int)*self.unknown_streamline # assign unknown label as default
        for k in range(n_streamlines):
            num_points = N.fromfile(f, dtype='<i4', count=1)[0]
            xyz_scalar = N.fromfile(f, dtype='<f4', count=num_points*(3+n_scalars)).reshape(num_points, 3+n_scalars)
            properties = N.fromfile(f, dtype='<f4', count=self.header['n_properties'][0])
            self.streamline.append(xyz_scalar)
            self.properties.append(properties)
            self.progress_meter(k, n_streamlines, 'Reading streamlines...')
            pass
        return
        
    def writeStreamlinesTrk(self, f):
        """Write streamlines to file in .trk format. Assumption: header has
        already been written.
        """
        n_scalars = self.header['n_scalars'][0]
        n_streamlines = self.header['n_count'][0]
        k = 0
        for streamline_id in self.streamline_id:
            num_points = N.array((self.getStreamline(streamline_id)).shape[0], dtype='<i4')
            num_points.tofile(f)
            xyz_scalar = N.array(self.getStreamline(streamline_id), dtype='<f4')
            xyz_scalar.tofile(f)
            properties = N.array(self.getProperties(streamline_id), dtype='<f4')
            properties.tofile(f)
            self.progress_meter(k, n_streamlines, 'Writing streamlines...')
            k += 1
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
        k = 0
        for streamline_id in self.streamline_id:
            xyz = self.getStreamline(streamline_id)
            ijk = self.mm2voxel(xyz)
            for i in range(xyz.shape[0]):
                try:
                    self.voxel2streamlines[tuple(ijk[i,:])].append(streamline_id)
                except KeyError:
                    self.voxel2streamlines[tuple(ijk[i,:])] = [streamline_id]
                    pass
                pass
            self.progress_meter(k, n_streamlines, 'Mapping voxels to streamlines...')
            k += 1
            pass
        n_voxels = len(self.voxel2streamlines.keys())
        # Now transform each list of IDs in an array of IDs:
        for n, ijk in enumerate(self.voxel2streamlines.keys()):
            self.voxel2streamlines[ijk] = N.array(self.voxel2streamlines[ijk])
            self.progress_meter(n, n_voxels, 'Converting lists to arrays...')
            pass
        return self.voxel2streamlines

    def fastCopy(self):
        """Create a copy of the current object avoiding the self.fiber
        list, which could be quite expensive.
        """
        new_streamlines = Streamlines()
        new_streamlines.header = copy.deepcopy(self.header)
        new_streamlines.num_label = copy.deepcopy(self.num_label)
        new_streamlines.label_num = copy.deepcopy(self.label_num)
        return new_streamlines

    def getStreamline(self, streamline_id):
        return self.streamline[self.streamlineid_pos[streamline_id]]

    def getProperties(self, streamline_id):
        return self.properties[self.streamlineid_pos[streamline_id]]

    def getStreamlineLabel(self, streamline_id):
        return self.streamline_label[self.streamlineid_pos[streamline_id]]

    def selectStreamlines(self, streamlines_ids):
        """Given a list of streamlines ID returns a new Streamline
        object with just those streamlines and an header equivalent to
        the current one except for header['n_count'], which is
        adjusted to the actual number of fibers.
        """
        new_streamlines = self.fastCopy()
        new_streamlines.streamline = [self.getStreamline(streamline_id) for streamline_id in streamlines_ids]
        new_streamlines.properties = [self.getProperties(streamline_id) for streamline_id in streamlines_ids]
        new_streamlines.streamline_id = list(streamlines_ids)
        new_streamlines.streamlineid_pos = dict(zip(streamlines_ids, range(len(streamlines_ids)))) # this remaps streamline_id to the correct new positions
        new_streamlines.streamline_label = N.array([self.getStreamlineLabel(streamline_id) for streamline_id in streamlines_ids])
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
        # streamlines_ids = N.unique(N.hstack([self.voxel2streamlines[tuple(v)] for v in voxels]))
        # ALTERNATIVE:
        # streamlines_ids = N.unique(N.hstack([self.voxel2streamlines[i,j,k] for i,j,k in voxels]))
        tmp = []
        for i,j,k in voxels:
            try:
                tmp.append(self.voxel2streamlines[i,j,k])
            except KeyError: # some voxels are empty...
                pass
            pass
        if tmp!=[]:
            streamlines_ids = N.unique(N.hstack(tmp))
            return self.selectStreamlines(streamlines_ids)
        return None

    def selectStreamlinesFromSlice(self, x=None, y=None, z=None):
        """Select streamlines specifying a slice (e.g., y=5)
        or a combination of slices (e.g., x=2, y=5).
        """
        # Build the list of voxels coordinates related to the desired slice:
        xx, yy, zz = 1, 1, 1
        if x==None: xx = self.header['dim'][0]; x = 0
        if y==None: yy = self.header['dim'][1]; y = 0
        if z==None: zz = self.header['dim'][2]; z = 0
        voxels = N.indices((xx,yy,zz)).reshape(3,xx*yy*zz).T + N.array([x,y,z])
        return self.selectStreamlinesFromVoxels(voxels)

    def getVolume(self, count=False):
        """Return a volume where a voxel is 1 if at least one fiber
        crosses it, otherwise 0.
        """
        volume = N.zeros(self.header['dim'], dtype='i')
        if (self.voxel2streamlines is not None) and count==False: # fast:
            ijk = N.array(self.voxel2streamlines.keys())
            volume[ijk[:,0],ijk[:,1],ijk[:,2]] = 1.0
            pass
        else: # slow but does not require voxel2streamlines:
            for xyz in self.streamline:
                ijk = self.mm2voxel(xyz)
                volume[ijk[:,0],ijk[:,1],ijk[:,2]] += 1
                pass
            pass
        if count:
            return volume
        return (volume>0).astype('i')

    def save(self, filename=None):
        """Dump Streamlines object to file. Assign a default filename if necessary.
        """
        if filename is None:
            self.dumpfilename = "Streamlines_"+self.filename+".pickle"
            pass
        print "Saving current instance to file",filename,"...",
        f = gzip.open(filename,'w')
        cPickle.dump(self, f, protocol=2)
        f.close()
        print os.stat(filename)[6],"bytes...",
        print "Done."
        return

    def setLabels(self, labeled_streamline_ids, labeled_streamline_numlabels, label_num=None, num_label=None):
        """Add labels to some streamlines
        """
        for i in range(len(labeled_streamline_ids)):
            lsid = labeled_streamline_ids[i]
            self.streamline_label[self.streamlineid_pos[lsid]] = labeled_streamline_numlabels[i]
            pass
        self.label_num = copy.deepcopy(label_num)
        self.num_label = copy.deepcopy(num_label)
        return
    
    
    def setNewVoxelSize(self, new_voxel_size = N.array([4.0, 4.0, 4.0])):
        """Set a new voxel size.
        """
        self.header['dim'] = (self.header['dim']*self.header['voxel_size']/new_voxel_size).astype('<h')
        print "New header['dim'] =", self.header['dim']
        self.header['voxel_size'] = new_voxel_size.astype('<f4')
        return

    def getVoxelSize(self):
        size = self.header['voxel_size']
        return size
 
    def buildStreamlineStreamlinesDict(self):
        """Build a dictionary mapping an input_streamline_id to an array of
        streamlines_ids that cross at least one of the voxels (of the grid)
        crossed by the input streamline. This means that the resluting streamlines
        are candidates to be clustered with the input streamline.

        Note that the size of the voxels can be changed at will with
        setNewVoxelSize().
        """
        n_streamlines = len(self.streamline)
        self.streamline2streamlines = {}
        counter = 0
        for streamline_id in self.streamline_id:
            xyz = self.getStreamline(streamline_id)
            ijk = self.mm2voxel(xyz)
            streamline_id_list = N.unique(N.hstack([self.voxel2streamlines[i,j,k] for i,j,k in ijk]))
            self.streamline2streamlines[streamline_id] = streamline_id_list
            self.progress_meter(counter, n_streamlines, "Building streamline to streamlines map...")
            counter += 1
            pass
        return self.streamline2streamlines


if __name__=="__main__":
    
    pass
