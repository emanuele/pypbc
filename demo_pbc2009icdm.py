# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

"""Demo program to extract main tracts from PBC2009ICDM data.

Copyright (c) 2009 Emanuele Olivetti <emanuele_AT_relativita.com>

This library is free software; you can redistribute it and/or modify
it either under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.
"""

import numpy as np
import streamlines


def build_filename(brain, scan):
    """Build full paths and filenames to reach PBC2009ICDM files.
    """
    datadir = 'pbc2009icdm/brain'+str(brain)+'/'
    base = 'brain'+str(brain)+'_scan'+str(scan)
    filename_trk = base+'_fiber_track_mni.trk'
    filename_label = base+'_fiber_labels.txt'
    filename_id_label = base+'_bundle_ids.txt'
    return {'trk':datadir+filename_trk, 'label':datadir+filename_label,
            'id_label':datadir+filename_id_label}


def parse_labels(brain, scan):
    """Parse labels from PBC2009ICDM labels file and
    return it as a matrix.
    """
    a=np.fromfile(build_filename(brain, scan)['label'], dtype=int, sep='\t')
    return a.reshape(a.shape[0]/2,2)


def parse_ID_label(brain, scan):
    """Parse ID->label file and build and return a dictionary.
    """
    ID_label = {}
    for line in file(build_filename(brain, scan)['id_label']):
        ID,name = line.split('       ') # there are excetly 7 spaces between ID and name
        ID_label[int(ID)] = name.strip()
        pass
    return ID_label


if __name__=="__main__":

    brain = 1
    scan = 1

    # read trackvis file:
    s = streamlines.Streamlines()
    s.loadTrk(build_filename(brain,scan)['trk'])
    s.printHeaderTrk()

    # parse labels and IDs:
    label = parse_labels(brain,scan)
    label_name = parse_ID_label(brain,scan)

    print "Extracting main tracts:"

    # create a trackvis file for each main tract:
    label_streamlineIDs = {}
    for lab in label_name.keys():
        label_streamlineIDs[lab] = label[(label[:,1]==lab),0]
        print lab, ')', label_name[lab],':' , label_streamlineIDs[lab].size, 'streamlines'
        ss = s.selectStreamlines(label_streamlineIDs[lab])
        ss.printHeaderTrk()
        tract_filename = 'brain'+str(brain)+'_scan'+str(scan)+'_fiber_track_mni_'+label_name[lab].replace(' ','_')+'.trk'
        ss.saveTrk(tract_filename)
        pass
    
