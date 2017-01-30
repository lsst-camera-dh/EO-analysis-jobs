#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
from __future__ import print_function
import os
import glob
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

for sensor_id in raft.sensor_names:
    fe55_files = siteUtils.dependency_glob('S*/%s_fe55_fe55_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                           description='Fe55 files:')
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                           description='Bias files:')
    nf = len(bias_files)
    mean_bias_file = '%(sensor_id)s_mean_bias_%(nf)i.fits' % locals()
    imutils.fits_mean_file(bias_files, mean_bias_file)

    #
    # Create a png zoom of the upper right corner of segment 1 for an Fe55
    # exposure for inclusion in the test report.
    #
    print("processing fe55_zoom:", fe55_files[0])
    sensorTest.fe55_zoom(fe55_files[0], size=250, amp=1)
    plt.savefig('%(sensor_id)s_fe55_zoom.png' % locals())

    #
    # Perform analysis of 9-pixel statistics for Fe55 charge clusters.
    #
    try:
        pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=sensor_id)
        pixel_stats.pixel_hists(pix0='p3', pix1='p5')
        plt.savefig('%(sensor_id)s_fe55_p3_p5_hists.png' % locals())

        pixel_stats.pixel_diff_profile(pixel_coord='x', pix0='p3', pix1='p5')
        plt.savefig('%(sensor_id)s_fe55_p3_p5_profiles.png' % locals())

        pixel_stats.apflux_profile()
        plt.savefig('%(sensor_id)s_fe55_apflux_serial.png' % locals())

        pixel_stats.apflux_profile(pixel_coord='y')
        plt.savefig('%(sensor_id)s_fe55_apflux_parallel.png' % locals())
    except Exception as eobj:
        print("Exception raised while creating pixel statistics plots:")
        print(str(eobj))
        print("Skipping these plots.")

    # Roll-off defects mask needs an input file to get the vendor
    # geometry, and will be used for all analyses.
    rolloff_mask_file = '%s_rolloff_defects_mask.fits' % sensor_id
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    task = sensorTest.Fe55Task()
    task.config.temp_set_point = -100.
    task.run(sensor_id, fe55_files, (rolloff_mask_file,), accuracy_req=0.01)

    # Fe55 gain and psf analysis results plots for the test report.
    results_file = '%s_eotest_results.fits' % sensor_id
    fe55_file = glob.glob('%s_psf_results*.fits' % sensor_id)[0]
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    try:
        png_file = '%s_fe55_dists.png' % sensor_id
        plots.fe55_dists(fe55_file=fe55_file)
        plt.savefig(png_file)
    except Exception as eobj:
        print("Exception raised while creating %s:" % png_file)
        print(eobj)

    try:
        png_file = '%s_psf_dists.png' % sensor_id
        plots.psf_dists(fe55_file=fe55_file)
        plt.savefig(png_file)
    except Exception as eobj:
        print("Exception raised while creating %s:" % png_file)
        print(eobj)

    try:
        png_file = '%s_gains.png' % sensor_id
        plots.gains()
        plt.savefig(png_file)
    except Exception as eobj:
        print("Exception raised while creating %s:" % png_file)
        print(eobj)

    sensorTest.plot_flat(mean_bias_file,
                         title='%s, mean bias frame' % sensor_id)
    plt.savefig('%s_mean_bias.png' % sensor_id)
