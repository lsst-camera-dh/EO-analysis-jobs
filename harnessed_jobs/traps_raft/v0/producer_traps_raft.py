#!/usr/bin/env python
"""
Producer script for raft-level traps analysis.
"""
from __future__ import print_function
import os
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from multiprocessor_execution import sensor_analyses

def run_trap_task(sensor_id):
    """Run the traps analysis for a single sensor."""
    trap_file = siteUtils.dependency_glob('S*/%s_trap_ppump_*.fits' % sensor_id,
                                          jobname=siteUtils.getProcessName('ppump_raft_acq'),
                                          description='Trap file:')[0]
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    # Omit rolloff defects mask since a trap in the rolloff edge region can
    # affect the entire column.
    mask_files = [item for item in mask_files
                  if item.find('rolloff_defects') == -1]
    print("Using mask files:")
    for mask_file in mask_files:
        print("  " + mask_file)

    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.TrapTask()
    task.run(sensor_id, trap_file, mask_files, gains)

if __name__ == '__main__':
    if os.environ.get('LCATR_SKIP_TRAPS_ANALYSIS', 'False') == 'True':
        pass
    else:
        sensor_analyses(run_trap_task)
