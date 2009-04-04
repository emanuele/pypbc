"""Plotting streamlines with mayavi2.

Note that you need a recent mayavi2 in order to run this. If you
use Ubuntu Hardy Heron 8.04 LTS like me, have a look here:
http://gael-varoquaux.info/blog/?p=70
and use these repositories:
https://launchpad.net/~gael-varoquaux/+archive
Note that mayavi2 shipped with more recent releases of Ubuntu works
fine.


Copyright (c) 2009 Emanuele Olivetti <emanuele_AT_relativita.com>

This library is free software; you can redistribute it and/or modify
it either under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.
"""

from enthought.mayavi import mlab
import streamlines


if __name__=="__main__":

    # filename = "DeterministicTractography/QBALLRecon/hardiO10.trk"
    filename = "hardiO10.trk_cross_streamline_id_20.trk"

    
    try:
        s
    except:
        s = streamlines.Streamlines()
        print
        print "file:", filename
        s.loadTrk(filename)
        s.printHeaderTrk()
        pass

    for stream, tmp in s.streamline:
        mlab.plot3d(stream[:,0], stream[:,1], stream[:,2], tube_radius=0.2)
        pass

    
    mlab.show()
