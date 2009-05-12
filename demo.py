# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""Demo program playing with fibers/streamlines.

Copyright (c) 2009 Emanuele Olivetti <emanuele_AT_relativita.com>

This library is free software; you can redistribute it and/or modify
it either under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.
"""

import streamlines
import sys
import numpy as N

if __name__=="__main__":
    
    print "This simple program reads a TrackVis .trk file, parse it, build"
    print "structures to represent streamlines as Python list of arrays"
    print "and then saves structures in TrackVis .trk file format."
    print "The resulting file is expected to be identical to the original."
    print "As a further step a dictionary, mapping voxel to streamlines, is built"
    print "and some examples using it are shown."
    
    
    s = streamlines.Streamlines()

    # filename = "dsi.trk"
    # filename = "dti.trk"
    # filename = "hardiO10.trk"
    filename = "hardiO10.trk_cross_streamline_id_2000.trk"

    if sys.argv[-1].endswith(".trk"): filename = sys.argv[-1]
    
    print
    print "file:", filename
    s.loadTrk(filename)
    s.printHeaderTrk()
    
    print
    streamline_id = 10
    print "Example: streamline_id=",streamline_id
    print s.get_streamline(streamline_id)
    print "Convert points from mm to voxel coordinates:"
    Vxyz = s.mm2voxel(s.get_streamline(streamline_id))
    print Vxyz
    print "Convert back and check whether differences are less than grid size...",
    assert(((s.voxel2mm(Vxyz) - s.get_streamline(streamline_id))<s.header['voxel_size']).all())
    print "OK."
    
    print
    filename2 = filename+"_COPY.trk"
    print "Saving to:", filename2
    s.saveTrk(filename2)
    
    print
    print "Building voxel2streamlines dictionary:"
    s.buildVoxelStreamlinesDict()
    voxel = s.voxel2streamlines.keys()[0]
    print "Example: streamlines crossing voxel", voxel
    print s.voxel2streamlines[voxel]
    
    print
    x = s.header['dim'][0] / 2
    print "Example: counting streamlines crossing plane x =", x
    counter = 0
    for y in range(s.header['dim'][1]):
        for z in range(s.header['dim'][2]):
            try:
                counter += s.voxel2streamlines[(x,y,z)].size
            except KeyError:
                pass
            pass
        pass
    print "Number of streamlines:", counter
    
    print
    streamline_id = 20
    print "Which streamlines cross (the voxels of) streamline with ID="+str(streamline_id)+" ?"
    xyz = s.get_streamline(streamline_id)
    ijk = s.mm2voxel(xyz)
    s2 = s.selectStreamlinesFromVoxels(ijk)
    print len(s2.streamline), "streamlines."

    # WARNING: s2 is a subset of the fibers of s but we lost the
    # streamline_id information since it is not stored anywhere!!! It
    # would be better to store IDs somewhere...
    
    print
    print "Saving .trk file with just the previous list of streamlines."
    filename3 = filename+'_cross_streamline_id_'+str(streamline_id)+'.trk'
    print "Saving to:", filename3
    s2.saveTrk(filename3)

    volume2 = s2.getVolume()
    print volume2.sum(), "voxels crossed by those streamlines."

    print
    print "Building voxel2streamlines dictionary:"
    s2.buildVoxelStreamlinesDict()
    volume2 = s2.getVolume()
    print volume2.sum(), "voxels crossed by those streamlines."
    
    volume3 = s2.getVolume(count=True)
    print volume3.max(), "is the max number of streamlines crossing a single voxel."
    print "this voxel is", N.unravel_index(volume3.argmax(), volume3.shape)

    print
    x = 50
    print "Selecting streamlines crossing slice x =",x
    s3 = s.selectStreamlinesFromSlice(x=x)
    print len(s3.streamline), "streamlines."

    # WARNING: I'm trying to select the same slice in TrackVis
    # using hardiO10.trk but I get 41065 streamlines there
    # (Skip=0%) of the instead of the 35571 I get here!!
    # Note that I compare slice x=50 here and slice x=51 in
    # TrackVis, according to the different numbering scheme.
