#!/usr/bin/env python
"""
Validator script for raft-level dark current analysis.
"""
import lsst.eotest.sensor as sensorTest
import lcatr.schema
import siteUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()
    results_file = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(results_file)

    amps = data['AMP']
    dc95s = data['DARK_CURRENT_95']
    dcmeds = data['DARK_CURRENT_MEDIAN']
    for amp, dc95, dcmed in zip(amps, dc95s, dcmeds):
        results.append(lcatr.schema.valid(lcatr.schema.get('dark_current_raft'),
                                          amp=amp, dark_current_95CL=dc95,
                                          dark_current_median=dcmed,
                                          slot=slot, sensor_id=sensor_id))

    # Persist the png files.
    metadata = dict(CCD_MANU=ccd_vendor, LSST_NUM=sensor_id,
                    TESTTYPE='DARK', TEST_CATEGORY='EO')
    results.extend(siteUtils.persist_png_files('%s*.png' % sensor_id,
                                               sensor_id, folder=slot,
                                               metadata=metadata))

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
