#!/usr/bin/env python
"""
Script to collect the EO analysis results for each sensor and write
an eotest_report.fits file.  Also create raft-level mosaics.
"""
from __future__ import print_function
import lsst.eotest.sensor as sensorTest
from lcatr.harness.helpers import dependency_glob
import eotestUtils
import siteUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

for sensor_id in raft.sensor_names:
    results_file = '%s_eotest_results.fits' % sensor_id

    # Use the mean bias file to determine the maximum number of
    # active pixels for the image quality statistics.
    bias_file = siteUtils.dependency_glob('S*/%s_*mean_bias*.fits' % sensor_id,
                                          jobname='fe55_raft_analysis',
                                          description="Mean bias files:")[0]
    total_num, rolloff_mask = sensorTest.pixel_counts(bias_file)

    # Aggregate information from summary.lims files into a final
    # EOTestResults output file.
    repackager = eotestUtils.JsonRepackager()
    repackager.eotest_results.add_ccd_result('TOTAL_NUM_PIXELS', total_num)
    repackager.eotest_results.add_ccd_result('ROLLOFF_MASK_PIXELS',
                                             rolloff_mask)
    summary_files = dependency_glob('summary.lims')
    repackager.process_files(summary_files)
    repackager.write(results_file)
