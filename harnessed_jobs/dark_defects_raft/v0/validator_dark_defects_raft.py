#!/usr/bin/env python
"""
Validator script for raft-level dark defects analysis.
"""
import glob
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

    mask_file = '%s_dark_pixel_mask.fits' % sensor_id
    eotestUtils.addHeaderData(mask_file, LSST_NUM=sensor_id,
                              TESTTYPE='SFLAT_500',
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor)
    results.append(siteUtils.make_fileref(mask_file))

    bias_frame = glob.glob(f'{sensor_id}_*_median_bias.fits')[0]
    eotestUtils.addHeaderData(bias_frame, DATE=eotestUtils.utc_now_isoformat())
    results.append(siteUtils.make_fileref(bias_frame, folder=slot))

    superflat = '%s_median_sflat.fits' % sensor_id
    eotestUtils.addHeaderData(superflat, DATE=eotestUtils.utc_now_isoformat())
    results.append(siteUtils.make_fileref(superflat))

    eotest_results = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(eotest_results)
    amps = data['AMP']
    npixels = data['NUM_DARK_PIXELS']
    ncolumns = data['NUM_DARK_COLUMNS']
    for amp, npix, ncol in zip(amps, npixels, ncolumns):
        results.append(lcatr.schema.valid(lcatr.schema.get('dark_defects_raft'),
                                          amp=amp,
                                          dark_pixels=npix,
                                          dark_columns=ncol,
                                          slot=slot,
                                          sensor_id=sensor_id))
    # Persist the png files.
    metadata = dict(CCD_MANU=ccd_vendor, LSST_NUM=sensor_id,
                    TESTTYPE='SFLAT_500', TEST_CATEGORY='EO')
    results.extend(siteUtils.persist_png_files('%s*.png' % sensor_id,
                                               sensor_id, folder=slot,
                                               metadata=metadata))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
