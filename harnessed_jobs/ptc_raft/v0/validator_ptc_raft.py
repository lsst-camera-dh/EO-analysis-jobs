#!/usr/bin/env ipython
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
    columns = (data['AMP'], data['PTC_GAIN'], data['PTC_GAIN_ERROR'],
               data['PTC_A00'], data['PTC_A00_ERROR'],
               data['PTC_NOISE'], data['PTC_NOISE_ERROR'],
               data['PTC_TURNOFF'])
    schema = lcatr.schema.get('ptc_raft')
    for (amp, ptc_gain, ptc_gain_error, ptc_a00, ptc_a00_error,
         ptc_noise, ptc_noise_error, ptc_turnoff) in zip(*columns):
        results.append(lcatr.schema.valid(
            schema, amp=amp,
            ptc_gain=ptc_gain, ptc_gain_error=ptc_gain_error,
            ptc_a00=ptc_a00, ptc_a00_error=ptc_a00_error,
            ptc_noise=ptc_noise, ptc_noise_error=ptc_noise_error,
            ptc_turnoff=ptc_turnoff, slot=slot, sensor_id=sensor_id))

    # Persist the png files.
    metadata = dict(CCD_MANU=ccd_vendor, LSST_NUM=sensor_id,
                    TESTTYPE='FLAT', TEST_CATEGORY='EO')
    results.extend(siteUtils.persist_png_files('%s*.png' % sensor_id,
                                               sensor_id, folder=slot,
                                               metadata=metadata))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
