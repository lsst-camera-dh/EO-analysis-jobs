#!/usr/bin/env python
"""
Producer script for raft-level CTE analysis.
"""
from __future__ import print_function
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
import camera_components

siteUtils.aggregate_job_ids()

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

for sensor_id in raft.sensor_names:
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)
    # Omit rolloff defects mask since it would mask some of the edges used
    # in the eper method.
    mask_files = [item for item in mask_files if
                  item.find('rolloff_defects') == -1]
    print("Using mask files:")
    for mask_file in mask_files:
        print("  " + mask_file)

    sflat_high_files = \
        siteUtils.dependency_glob('S*/%s_sflat_500_flat_H*.fits' % sensor_id,
                                  jobname=siteUtils.getProcessName('sflat_raft_acq'),
                                  description='Superflat high flux files:')

    task = sensorTest.CteTask()
    task.run(sensor_id, sflat_high_files, flux_level='high', gains=gains,
             mask_files=mask_files)

    sflat_low_files = \
        siteUtils.dependency_glob('S*/%s_sflat_500_flat_L*.fits' % sensor_id,
                                  jobname=siteUtils.getProcessName('sflat_raft_acq'),
                                  description='Superflat low flux files:')
    task.run(sensor_id, sflat_low_files, flux_level='low', gains=gains,
             mask_files=mask_files)
