#!/usr/bin/env python
"""
Producer script for raft-level PRNU analysis.
"""
from __future__ import print_function
import sys
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
import camera_components

siteUtils.aggregate_job_ids()

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

for sensor_id in raft.sensor_names:
    lambda_files = siteUtils.dependency_glob('S*/%s_lambda_flat_*.fits' % sensor_id,
                                             jobname=siteUtils.getProcessName('qe_raft_acq'),
                                             description='Lambda files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    correction_image = eotestUtils.getIlluminationNonUniformityImage()
    if correction_image is None:
        print()
        print("WARNING: The correction image file is not given in")
        print("config/%s/eotest_calibrations.cfg." % siteUtils.getSiteName())
        print("No correction for non-uniform illumination will be applied.")
        print()
        sys.stdout.flush()

    task = sensorTest.PrnuTask()
    task.run(sensor_id, lambda_files, mask_files, gains, correction_image)

