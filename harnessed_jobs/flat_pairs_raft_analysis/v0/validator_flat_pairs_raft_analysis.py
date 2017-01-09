#!/usr/bin/env python
"""
Validator script for raft-level flat pairs analysis.
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

    det_resp_data = '%s_det_response.fits' % sensor_id
    eotestUtils.addHeaderData(det_resp_data, LSST_NUM=sensor_id,
                              TESTTYPE='FLAT',
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor)
    results.append(siteUtils.make_fileref(det_resp_data, folder=slot))

    results_file = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(results_file)
    amps = data['AMP']
    full_well_data = data['FULL_WELL']
    max_frac_dev_data = data['MAX_FRAC_DEV']

    for amp, full_well, max_frac_dev in zip(amps, full_well_data,
                                            max_frac_dev_data):
        results.append(lcatr.schema.valid(lcatr.schema.get('flat_pairs_raft'),
                                          amp=amp, full_well=full_well,
                                          max_frac_dev=max_frac_dev,
                                          slot=slot,
                                          sensor_id=sensor_id))

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
