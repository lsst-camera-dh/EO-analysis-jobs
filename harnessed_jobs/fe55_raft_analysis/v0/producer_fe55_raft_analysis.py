#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
from __future__ import print_function
import os
# This is needed so that pyplot can write to .matplotlib
os.environ['MPLCONFIGDIR'] = os.curdir
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import siteUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

print("fe55_raft_acq process name", siteUtils.getProcessName('fe55_raft_acq'))

for sensor_id in raft.sensor_names:
    fe55_files = siteUtils.dependency_glob('S*/%s_fe55_fe55_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                           description='Fe55 files:')[:1]
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                           description='Bias files:')
    nf = len(bias_files)
    outfile = '%(sensor_id)s_mean_bias_%(nf)i.fits' % locals()
    imutils.fits_mean_file(bias_files, outfile)

    #
    # Create a png zoom of the upper right corner of segment 1 for an Fe55
    # exposure for inclusion in the test report
    #
    print("processing fe55_zoom:", fe55_files[0])
    sensorTest.fe55_zoom(fe55_files[0], size=250, amp=1)
    plt.savefig('%(sensor_id)s_fe55_zoom.png' % locals())

    #
    # Perform analysis of 9-pixel statistics for Fe55 charge clusters.
    #
    pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=sensor_id)
    pixel_stats.pixel_hists(pix0='p3', pix1='p5')
    plt.savefig('%(sensor_id)s_fe55_p3_p5_hists.png' % locals())

    pixel_stats.pixel_diff_profile(pixel_coord='x', pix0='p3', pix1='p5')
    plt.savefig('%(sensor_id)s_fe55_p3_p5_profiles.png' % locals())

    pixel_stats.apflux_profile()
    plt.savefig('%(sensor_id)s_fe55_apflux_serial.png' % locals())

    pixel_stats.apflux_profile(pixel_coord='y')
    plt.savefig('%(sensor_id)s_fe55_apflux_parallel.png' % locals())

    # Roll-off defects mask needs an input file to get the vendor
    # geometry, and will be used for all analyses.
    rolloff_mask_file = '%s_rolloff_defects_mask.fits' % sensor_id
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    task = sensorTest.Fe55Task()
    task.run(sensor_id, fe55_files, (rolloff_mask_file,), accuracy_req=0.01)
