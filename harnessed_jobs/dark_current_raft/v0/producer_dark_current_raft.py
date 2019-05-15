#!/usr/bin/env python
"""
Producer script for raft-level dark current analysis.
"""
from __future__ import absolute_import
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from multiprocessor_execution import sensor_analyses

def run_dark_current_task(sensor_id):
    "Single sensor execution of dark current analysis."
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    dark_files = siteUtils.dependency_glob('S*/%s_dark_dark_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('dark_raft_acq'),
                                           description='Dark files:')
    bias_frame = siteUtils.dependency_glob('%s_mean_bias*.fits' % sensor_id,
                                           description='Super bias frame:')[0]
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.DarkCurrentTask()
    task.config.temp_set_point = -100.
    dark_curr_pixels, dark95s \
        = task.run(sensor_id, dark_files, mask_files, gains,
                   bias_frame=bias_frame)

    results_file \
        = siteUtils.dependency_glob('%s_eotest_results.fits' % sensor_id,
                                    jobname='read_noise_raft')[0]
#    eo_results = sensorTest.EOTestResults(results_file)
#    read_noise = dict(pair for pair in zip(eo_results['AMP'],
#                                           eo_results['TOTAL_NOISE']))
#
#    siteUtils.make_png_file(sensorTest.total_noise_histograms,
#                            '%s_total_noise_hists.png' % file_prefix,
#                            dark_curr_pixels, read_noise, dark95s,
#                            exptime=16, title=sensor_id)
#
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)
    siteUtils.make_png_file(plots.total_noise, '%s_noise.png' % file_prefix,
                            dark95s=dark95s)

if __name__ == '__main__':
    sensor_analyses(run_dark_current_task)
