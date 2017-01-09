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

    ptc_results = '%s_ptc.fits' % sensor_id
    eotestUtils.addHeaderData(ptc_results, LSST_NUM=sensor_id, TESTTYPE='FLAT',
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor)

    results.append(siteUtils.make_fileref(ptc_results, folder=slot))

    results_file = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(results_file)
    amps = data['AMP']
    ptc_gains = data['PTC_GAIN']
    ptc_gain_errors = data['PTC_GAIN_ERROR']
    for amp, gain, gain_error in zip(amps, ptc_gains, ptc_gain_errors):
        results.append(lcatr.schema.valid(lcatr.schema.get('ptc_raft'),
                                          amp=amp, ptc_gain=gain,
                                          ptc_gain_error=gain_error,
                                          slot=slot,
                                          sensor_id=sensor_id))

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
