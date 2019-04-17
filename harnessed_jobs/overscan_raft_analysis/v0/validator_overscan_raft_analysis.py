#!/usr/bin/env python
"""
Validator script for raft-level PTC analysis.
"""
import lsst.eotest.sensor as sensorTest
import lcatr.schema
import siteUtils
import eotestUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()

    overscan_results = '%s_overscan_results.fits' % sensor_id
    eotestUtils.addHeaderData(overscan_results, LSST_NUM=sensor_id, TESTTYPE='FLAT',
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor)

    results.append(siteUtils.make_fileref(overscan_results, folder=slot))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
