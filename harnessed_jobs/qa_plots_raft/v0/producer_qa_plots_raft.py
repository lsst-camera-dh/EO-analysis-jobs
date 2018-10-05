#!/usr/bin/env python
"""
QA plots for raft-level acquisitions.
"""
from __future__ import print_function
import os
import glob
from collections import OrderedDict
from lcatr.harness.helpers import dependency_glob
import siteUtils
from eo_acq_qa import RaftTrendingObjects
import camera_components

def dirname_dependencyGlob(ccd_vendor, **kwds):
    """
    Return the directory path with the FITS output files for the
    specified jobname.
    """
    if 'jobname' in kwds:
        kwds['jobname'] = siteUtils.getProcessName(kwds['jobname'])
    file0 = dependency_glob('S*/%(ccd_vendor)s*.fits' % locals(), **kwds)[0]
    # Apply os.path.dirname twice to omit both slot folder and file
    # basename.
    dirname = os.path.dirname(os.path.dirname(file0))
    if dirname == '':
        dirname = '.'
    return os.path.join(dirname, 'S*')

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

datasets = OrderedDict([('FE55', 'fe55_raft_acq'),
                        ('DARK', 'dark_raft_acq'),
                        ('FLAT', 'flat_pair_raft_acq'),
                        ('PPUMP', 'ppump_raft_acq'),
                        ('SFLAT', 'sflat_raft_acq'),
                        ('QE', 'qe_raft_acq')])

ccd_vendor = raft.sensor_names[0].split('-')[0]

QA_trender = RaftTrendingObjects(raft.sensor_names)
for test_type, jobname in datasets.items():
    dirname = dirname_dependencyGlob(ccd_vendor, jobname=jobname)
    print(dirname)
    QA_trender.processDirectory(dirname, test_type)

QA_trender.plot(raft_id)

file_prefix = '%s_%s' % (raft_id, siteUtils.getRunNumber())
for png_file in glob.glob('*QA*.png'):
    os.rename(png_file, png_file.replace(raft_id, file_prefix))
