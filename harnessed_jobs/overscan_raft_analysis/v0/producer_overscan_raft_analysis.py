#!/usr/bin/env python
"""
Producer script for raft-level overscan analysis.
"""
from __future__ import print_function
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from multiprocessor_execution import sensor_analyses

def run_overscan_task(sensor_id):

    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    flat_files = siteUtils.dependency_glob('S*/%s_flat*flat1_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('flat_pair_raft_acq'),
                                           description='Flat files:')
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.OverscanTask()
    task.run(sensor_id, flat_files, gains) # Needs bias file

    overscan_file = '%s_overscan_results.fits' % sensor_id
#    plots = sensorTest.EOTestPlots(sensor_id, overscan_file=overscan_file)

#    siteUtils.make_png_file(plots.overscan,
#                            '%s_linearity.png' % file_prefix)

if __name__ == '__main__':
    sensor_analyses(run_overscan_task)
