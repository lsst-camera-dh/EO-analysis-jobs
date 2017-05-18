#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
from __future__ import print_function
import glob
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import siteUtils
from multiprocessor_execution import sensor_analyses

def run_fe55_task(sensor_id):
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
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
    siteUtils.make_png_file(sensorTest.fe55_zoom,
                            '%(file_prefix)s_fe55_zoom.png' % locals(),
                            fe55_files[0], size=250, amp=1)

    #
    # Perform analysis of 9-pixel statistics for Fe55 charge clusters.
    #
    try:
        pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=sensor_id)

        siteUtils.make_png_file(pixel_stats.pixel_hists,
                                '%s_fe55_p3_p5_hists.png' % file_prefix,
                                pix0='p3', pix1='p5')

        siteUtils.make_png_file(pixel_stats.pixel_diff_profile,
                                '%s_fe55_p3_p5_profiles.png' % file_prefix,
                                pixel_coord='x', pix0='p3', pix1='p5')

        siteUtils.make_png_file(pixel_stats.apflux_profile,
                                '%s_fe55_apflux_serial.png' % file_prefix)

        siteUtils.make_png_file(pixel_stats.apflux_profile,
                                '%s_fe55_apflux_parallel.png' % file_prefix,
                                pixel_coord='y')

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
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    siteUtils.make_png_file(plots.gains,
                            '%s_gains.png' % file_prefix)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_mean_bias.png' % file_prefix,
                            mean_bias_file,
                            title='%s, mean bias frame' % sensor_id)

    fe55_file = glob.glob('%s_psf_results*.fits' % sensor_id)[0]
    siteUtils.make_png_file(plots.fe55_dists,
                            '%s_fe55_dists.png' % file_prefix,
                            fe55_file=fe55_file)

    siteUtils.make_png_file(plots.psf_dists,
                            '%s_psf_dists.png' % file_prefix,
                            fe55_file=fe55_file)

if __name__ == '__main__':
    sensor_analyses(run_fe55_task)
