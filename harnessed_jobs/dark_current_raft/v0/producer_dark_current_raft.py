#!/usr/bin/env python
"""
Producer script for raft-level dark current analysis.
"""
from __future__ import print_function
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from multiprocessor_execution import sensor_analyses

def run_dark_current_task(sensor_id):
    dark_files = siteUtils.dependency_glob('S*/%s_dark_dark_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('dark_raft_acq'),
                                           description='Dark files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.DarkCurrentTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, dark_files, mask_files, gains)

if __name__ == '__main__':
    sensor_analyses(run_dark_current_task)
