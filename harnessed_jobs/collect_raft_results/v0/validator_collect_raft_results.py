#!/usr/bin/env python
"""
Validator script for collect_raft_results job.
"""
from __future__ import print_function
import lcatr.schema
import siteUtils
import eotestUtils
import camera_components

results = []

run_number = siteUtils.getRunNumber()

md = siteUtils.DataCatalogMetadata(ORIGIN=siteUtils.getSiteName(),
                                   TEST_CATEGORY='EO',
                                   DATA_PRODUCT='EOTEST_RESULTS')

# Persist eotest_results files for each sensor.
raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()
    results_file = '%s_eotest_results.fits' % sensor_id
    eotestUtils.addHeaderData(results_file, LSST_NUM=sensor_id,
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor,
                              RUNNUM=run_number)
    results.append(lcatr.schema.fileref.make(results_file,
                                             metadata=md(CCD_MANU=ccd_vendor,
                                                         LSST_NUM=sensor_id,
                                                         SLOT=slot,
                                                         LsstId=raft_id)))

# Persist raft-level mosaics
# TBD

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
