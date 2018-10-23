#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
from __future__ import print_function
import glob
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import siteUtils
from camera_components import camera_info
from multiprocessor_execution import run_sensor_analysis_pool

def fe55_task(det_name):
    "Single sensor execution of the Fe55 analysis task."
    fe55_files \
        = siteUtils.dependency_glob('fe55_fe55_*/*_{}.fits'.format(det_name))
    bias_files \
        = siteUtils.dependency_glob('fe55_bias_*/*_{}.fits'.format(det_name))
    file_prefix = '{}_{}'.format(siteUtils.getRunNumber(), det_name)

    mean_bias_file = '{}_mean_bias.fits'.format(file_prefix)
    imutils.fits_mean_file(bias_files, mean_bias_file)

    pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=det_name)

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

    rolloff_mask_file = '%s_edge_rolloff_mask.fits' % file_prefix
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    task = sensorTest.Fe55Task()
    task.config.temp_set_point = -100.
    task.run(file_prefix, fe55_files, (rolloff_mask_file,), accuracy_req=0.01)

    # Fe55 gain and psf analysis results plots for the test report.
    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)

    siteUtils.make_png_file(plots.gains, '%s_gains.png' % file_prefix)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_mean_bias.png' % file_prefix,
                            mean_bias_file,
                            title='%s, mean bias frame' % det_name,
                            annotation='ADU/pixel, overscan-subtracted')

    fe55_file = glob.glob('%s_psf_results*.fits' % file_prefix)[0]
    siteUtils.make_png_file(plots.fe55_dists, '%s_fe55_dists.png' % file_prefix,
                            fe55_file=fe55_file)

    siteUtils.make_png_file(plots.psf_dists, '%s_psf_dists.png' % file_prefix,
                            fe55_file=fe55_file)

if __name__ == '__main__':
    det_names = camera_info.get_det_names()
    run_sensor_analysis_pool(fe55_task, det_names, processes=1)
