#!/usr/bin/env python
"""
Validator script for raft-level brighter-fatter analysis.
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

    bf_results = '%s_bf.fits' % sensor_id
    eotestUtils.addHeaderData(bf_results, LSST_NUM=sensor_id, TESTTYPE='FLAT',
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor)

    results.append(siteUtils.make_fileref(bf_results, folder=slot))

    results_file = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(results_file)
    amps = data['AMP']
    bf_xcorr = data['BF_XCORR']
    bf_xcorr_err = data['BF_XCORR_ERR']
    bf_ycorr = data['BF_YCORR']
    bf_ycorr_err = data['BF_YCORR_ERR']
    bf_mean = data['PTC_MEAN']
    columns = (data['AMP'], data['BF_XCORR'], data['BF_XCORR_ERR'],
               data['BF_YCORR'], data['BF_YCORR_ERR'], data['BF_MEAN'])
    schema = lcatr.schema.get('brighter_fatter_raft')
    for amp, bf_xcorr, bf_xcorr_err, bf_ycorr, bf_ycorr_err, bf_mean \
        in zip(*columns):
        results.append(lcatr.schema.valid(
            schema, amp=amp, bf_xcorr=bf_xcorr, bf_xcorr_err=bf_xcorr_err,
            bf_ycorr=bf_ycorr, bf_ycorr_err=bf_ycorr_err, bf_mean=bf_mean,
            slot=slot, sensor_id=sensor_id))

    # Persist the png files.
    metadata = dict(CCD_MANU=ccd_vendor, LSST_NUM=sensor_id,
                    TESTTYPE='FLAT', TEST_CATEGORY='EO')
    results.extend(siteUtils.persist_png_files('%s*.png' % sensor_id,
                                               sensor_id, folder=slot,
                                               metadata=metadata))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
