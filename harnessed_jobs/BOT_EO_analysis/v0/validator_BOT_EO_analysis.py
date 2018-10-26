#!/usr/bin/env python
"""
Validator script for BOT-level Fe55 analysis.
"""
from __future__ import print_function
import os
import glob
import numpy as np
import lcatr.schema
import siteUtils
import eotestUtils
import lsst.eotest.sensor as sensorTest
from camera_components import camera_info

def validate_fe55(results, det_names):
    """Validate and persist fe55 gain and psf results."""
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = '{}_{}'.format(siteUtils.getRunNumber(), det_name)

        # The output files from producer script.
        gain_file = '%(file_prefix)s_eotest_results.fits' % locals()
        psf_results_files \
            = glob.glob('%(file_prefix)s_psf_results*.fits' % locals())

        if not os.path.isfile(gain_file) or not psf_results_files:
            # Results for this detector are not available so skip it.
            continue
        psf_results = psf_results_files[0]

        rolloff_mask = '%(file_prefix)s_edge_rolloff_mask.fits' % locals()

        output_files = gain_file, psf_results, rolloff_mask

        # Add/update the metadata to the primary HDU of these files.
        for fitsfile in output_files:
            eotestUtils.addHeaderData(fitsfile, TESTTYPE='FE55',
                                      DATE=eotestUtils.utc_now_isoformat())
        results.extend([lcatr.schema.fileref.make(x) for x in output_files])

        # Persist the mean bias FITS file.
        bias_mean_file \
            = glob.glob('%(file_prefix)s_mean_bias.fits' % locals())[0]
        results.append(siteUtils.make_fileref(bias_mean_file))

        # Persist the png files.
        metadata = dict(TESTTYPE='FE55', TEST_CATEGORY='EO')
        results.extend(siteUtils.persist_png_files('%s*.png' % file_prefix,
                                                   det_name, metadata=metadata))

        data = sensorTest.EOTestResults(gain_file)
        amps = data['AMP']
        gain_data = data['GAIN']
        gain_errors = data['GAIN_ERROR']
        sigmas = data['PSF_SIGMA']
        for amp, gain_value, gain_error, sigma in zip(amps, gain_data,
                                                      gain_errors, sigmas):
            if not np.isfinite(gain_error):
                gain_error = -1
            results.append(lcatr.schema.valid(
                lcatr.schema.get('fe55_BOT_analysis'), amp=amp, gain=gain_value,
                gain_error=gain_error, psf_sigma=sigma, slot=slot, raft=raft))

    return results


def validate_read_noise(results, det_names):
    """Validate and persist read noise results."""
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = '{}_{}'.format(siteUtils.getRunNumber(), det_name)

        read_noise_file = '%s_eotest_results.fits' % file_prefix
        if not os.path.isfile(read_noise_file):
            print("read noise results not available for", det_name)
            continue
        data = sensorTest.EOTestResults(read_noise_file)
        amps = data['AMP']
        read_noise_data = data['READ_NOISE']
        system_noise_data = data['SYSTEM_NOISE']
        total_noise_data = data['TOTAL_NOISE']
        for amp, read_noise, system_noise, total_noise \
            in zip(amps, read_noise_data, system_noise_data, total_noise_data):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('read_noise_BOT'),
                amp=amp, read_noise=read_noise, system_noise=system_noise,
                total_noise=total_noise, slot=slot, raft=raft))

        files = glob.glob('%s_read_noise?*.fits' % file_prefix)
        for fitsfile in files:
            eotestUtils.addHeaderData(fitsfile, TESTTYPE='FE55',
                                      DATE=eotestUtils.utc_now_isoformat())

        data_products = [siteUtils.make_fileref(item) for item in files]
        results.extend(data_products)

        # Persist the png files.
        metadata = dict(TESTTYPE='FE55', TEST_CATEGORY='EO')
        results.extend(siteUtils.persist_png_files('%s*.png' % file_prefix,
                                                   det_name, metadata=metadata))

    # Persist the raft-level overscan correlation plots.
    metadata = dict(TESTTYPE='FE55', TEST_CATEGORY='EO')
    results.extend(siteUtils.persist_png_files('%s*.png' % file_prefix,
                                               raft, metadata=metadata))
    return results


def validate_bright_defects(results, det_names):
    """Validate and persist bright defects results."""
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = '{}_{}'.format(siteUtils.getRunNumber(), det_name)
        mask_file = '%s_bright_pixel_mask.fits' % file_prefix
        if not os.path.isfile(mask_file):
            print("mask file not found for {}. Skipping.".format(det_name))
            continue
        eotestUtils.addHeaderData(mask_file, TESTTYPE='DARK',
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(mask_file))

        medianed_dark = '%s_median_dark_bp.fits' % file_prefix
        eotestUtils.addHeaderData(medianed_dark,
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(medianed_dark))

        eotest_results = '%s_eotest_results.fits' % file_prefix
        data = sensorTest.EOTestResults(eotest_results)
        amps = data['AMP']
        npixels = data['NUM_BRIGHT_PIXELS']
        ncolumns = data['NUM_BRIGHT_COLUMNS']
        for amp, npix, ncol in zip(amps, npixels, ncolumns):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('bright_defects_BOT'),
                amp=amp, bright_pixels=npix, bright_columns=ncol,
                slot=slot, raft=raft))
        # Persist the png files.
        metadata = dict(TESTTYPE='DARK', TEST_CATEGORY='EO')
        results.extend(siteUtils.persist_png_files('%s*.png' % file_prefix,
                                                   det_name, metadata=metadata))
    return results


if __name__ == '__main__':
    det_names = camera_info.get_det_names()
    results = []
    results = validate_fe55(results, det_names)
    results = validate_read_noise(results, det_names)
    results = validate_bright_defects(results, det_names)

    results.extend(siteUtils.jobInfo())

    lcatr.schema.write_file(results)
    lcatr.schema.validate_file()
