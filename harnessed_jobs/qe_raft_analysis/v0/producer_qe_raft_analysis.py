#!/usr/bin/env python
"""
Producer script for raft-level QE analysis.
"""
from __future__ import print_function
import os
import sys
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

for sensor_id in raft.sensor_names:
    lambda_files = siteUtils.dependency_glob('S*/%s_lambda_flat_*.fits' % sensor_id,
                                             jobname=siteUtils.getProcessName('qe_raft_acq'),
                                             description='Lambda files:')

    pd_ratio_file = eotestUtils.getPhotodiodeRatioFile()
    if pd_ratio_file is None:
        message = ("The test-stand specific photodiode ratio file is " +
                   "not given in config/%s/eotest_calibrations.cfg."
                   % siteUtils.getSiteName())
        raise RuntimeError(message)

    correction_image = eotestUtils.getIlluminationNonUniformityImage()
    if correction_image is None:
        print()
        print("WARNING: The correction image file is not given in")
        print("config/%s/eotest_calibrations.cfg." % siteUtils.getSiteName())
        print("No correction for non-uniform illumination will be applied.")
        print()
        sys.stdout.flush()

    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.QeTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, lambda_files, pd_ratio_file, mask_files, gains,
             correction_image=correction_image)

    results_file \
        = siteUtils.dependency_glob('%s_eotest_results.fits' % sensor_id,
                                    jobname='fe55_raft_analysis',
                                    description='Fe55 results file')[0]
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    siteUtils.make_png_file(plots.qe,
                            '%s_qe.png' % sensor_id,
                            qe_file='%s_QE.fits' % sensor_id)

    wl_files \
        = siteUtils.dependency_glob('%s*_lambda_*.fits' % sensor_id,
                                    jobname=siteUtils.processName('qe_raft_acq'),
                                    description='Lambda files:')
    try:
        plots.flat_fields(os.path.dirname(wl_files[0]))
    except Exception as eobj:
        print("Exception raised while creating flat fields:")
        print(str(eobj))
