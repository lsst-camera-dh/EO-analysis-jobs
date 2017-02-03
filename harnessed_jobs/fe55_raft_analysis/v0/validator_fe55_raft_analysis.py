#!/usr/bin/env python
"""
Validator script for raft-level Fe55 analysis.
"""
from __future__ import print_function
import glob
import numpy as np
import lcatr.schema
import siteUtils
import eotestUtils
import lsst.eotest.sensor as sensorTest
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()
    # The output files from producer script.
    gain_file = '%(sensor_id)s_eotest_results.fits' % locals()
    psf_results = glob.glob('%(sensor_id)s_psf_results*.fits' % locals())[0]
    rolloff_mask = '%(sensor_id)s_rolloff_defects_mask.fits' % locals()

    output_files = gain_file, psf_results, rolloff_mask

    # Add/update the metadata to the primary HDU of these files.
    for fitsfile in output_files:
        eotestUtils.addHeaderData(fitsfile, LSST_NUM=sensor_id, TESTTYPE='FE55',
                                  DATE=eotestUtils.utc_now_isoformat(),
                                  CCD_MANU=ccd_vendor)

    #
    # Persist the mean bias FITS file.
    #
    bias_mean_file = glob.glob('%(sensor_id)s_mean_bias_*.fits' % locals())[0]
    results.append(siteUtils.make_fileref(bias_mean_file, folder=slot))

    # Persist the png files.
    metadata = dict(CCD_MANU=ccd_vendor, LSST_NUM=sensor_id,
                    TESTTYPE='FE55', TEST_CATEGORY='EO')
    results.extend(siteUtils.persist_png_files('%s*.png' % sensor_id,
                                               sensor_id, folder=slot,
                                               metadata=metadata))

    data = sensorTest.EOTestResults(gain_file)
    amps = data['AMP']
    gain_data = data['GAIN']
    gain_errors = data['GAIN_ERROR']
    sigmas = data['PSF_SIGMA']
    for amp, gain_value, gain_error, sigma in zip(amps, gain_data, gain_errors,
                                                  sigmas):
        if not np.isfinite(gain_error):
            gain_error = -1
        results.append(lcatr.schema.valid(lcatr.schema.get('fe55_raft_analysis'),
                                          amp=amp, gain=gain_value,
                                          gain_error=gain_error,
                                          psf_sigma=sigma,
                                          slot=slot,
                                          sensor_id=sensor_id))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
