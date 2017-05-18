#!/usr/bin/env python
"""
Producer script for raft-level dark defects analysis.
"""
from __future__ import print_function
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from multiprocessor_execution import sensor_analyses

def run_dark_pixels_task(sensor_id):
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    sflat_files = siteUtils.dependency_glob('S*/%s_sflat_500_flat_H*.fits' % sensor_id,
                                            jobname=siteUtils.getProcessName('sflat_raft_acq'),
                                            description='Superflat files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)

    task = sensorTest.DarkPixelsTask()
    task.run(sensor_id, sflat_files, mask_files)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_superflat_dark_defects.png' % file_prefix,
                            '%s_median_sflat.fits' % sensor_id,
                            title='%s, superflat for dark defects analysis' % sensor_id)

if __name__ == '__main__':
    sensor_analyses(run_dark_pixels_task)
