# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""Preliminary function to read and use brain regions from FreeSurfer
files.

Requires PyNIfTI:
http://niftilib.sourceforge.net/pynifti/

Copyright (c) 2009 Emanuele Olivetti <emanuele_AT_relativita.com>

This library is free software; you can redistribute it and/or modify
it either under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.
"""

import nifti
import numpy as N
import pylab as P

def parse_fs_lut(filename="FreeSurferColorLUT.txt"):
    """Parse FreeSurfer LookUp Table and return an equivalent
    dictionary.
    """
    f = open(filename)
    fs_lut = {}
    for line in f:
        token = line.split()
        if len(token)>=6:
            try:
                fs_lut[int(token[0])] = token[1]
                # print token[0], token[1]
            except ValueError:
                continue
            pass
        pass
    return fs_lut
                   

if __name__=="__main__":

    path = "PghBC2009/brain0/Structurals/"
    label_filename = "docs/FreeSurferColorLUT.txt"
    brain_filename = "DWISpace/brain0_asegonB0Anz.hdr"

    fs_lut = parse_fs_lut(path+label_filename)

    im = nifti.NiftiImage(path+brain_filename)
    volume = im.asarray()

    # Show one random slice:
    P.imshow(volume[44,:,:], interpolation='nearest')

    print "Regions labelled in this brain:"
    print "ID \tName \t\t#Voxels"
    for value in N.unique(volume):
        print value, "\t",fs_lut[value], "\t",(volume==value).sum()
        pass

    P.show()
