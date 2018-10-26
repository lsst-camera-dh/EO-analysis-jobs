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
    run = siteUtils.getRunNumber()
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = '{}_{}'.format(run, det_name)

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
        with open('fe55_task_png_files.txt', 'r') as input_:
            png_files = input_.readlines().strip().split()
        metadata = dict(TESTTYPE='FE55', TEST_CATEGORY='EO',
                        DETECTOR=det_name, RUN=run)
        results.extend(siteUtils.persist_png_files('', file_prefix,
                                                   png_files=png_files,
                                                   metadata=metadata))

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
    run = siteUtils.getRunNumber()
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = '{}_{}'.format(run, det_name)

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
        metadata = dict(DETECTOR=det_name, TESTTYPE='FE55', TEST_CATEGORY='EO',
                        RUN=run)
        filename = '%s_correlated_noise.png' % file_prefix
        results.extend(siteUtils.persist_png_files(filename, file_prefix,
                                                   metadata=metadata))

    # Persist the raft-level overscan correlation plots.
    for raft in camera_info.get_raft_names():
        metadata = dict(TESTTYPE='FE55', TEST_CATEGORY='EO', RAFT=raft, RUN=run)
        file_prefix = '{}_{}'.format(run, raft)
        filename = '%s_overscan_correlations.png' % file_prefix
        results.extend(siteUtils.persist_png_files(filename, file_prefix,
                                                   metadata=metadata))
    return results


def validate_bright_defects(results, det_names):
    """Validate and persist bright defects results."""
    run = siteUtils.getRunNumber()
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = '{}_{}'.format(run, det_name)
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
        # Persist the png file.
        metadata = dict(TESTTYPE='DARK', TEST_CATEGORY='EO', DETECTOR=det_name,
                        RUN=run)
        filename = '%s_medianed_dark.png' % file_prefix
        results.extend(siteUtils.persist_png_files(filename, file_prefix,
                                                   metadata=metadata))
    return results

def validate_dark_defects(results, det_names):
    """Validate and persist dark defects results."""
    run = siteUtils.getRunNumber()
    for det_name in det_names:
        file_prefix = '{}_{}'.format(run, det_name)
        mask_file = '%s_dark_pixel_mask.fits' % det_name
        eotestUtils.addHeaderData(mask_file, TESTTYPE='SFLAT_500',
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(mask_file))

        superflat = '%s_median_sflat.fits' % det_name
        eotestUtils.addHeaderData(superflat,
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(superflat))

        eotest_results = '%s_eotest_results.fits' % det_name
        data = sensorTest.EOTestResults(eotest_results)
        amps = data['AMP']
        npixels = data['NUM_DARK_PIXELS']
        ncolumns = data['NUM_DARK_COLUMNS']
        for amp, npix, ncol in zip(amps, npixels, ncolumns):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('dark_defects_BOT'),
                amp=amp, dark_pixels=npix, dark_columns=ncol,
                slot=slot, raft=raft))

        # Persist the png files.
        metadata = dict(DETECTOR=det_name, raft=raft,
                        TESTTYPE='SFLAT_500', TEST_CATEGORY='EO')
        filename = '%s_superflat_dark_defects.png' % file_prefix
        results.extend(siteUtils.persist_png_files(filename, file_prefix,
                                                   metadata=metadata))
    return results

def validate_traps(results, det_names):
    """Validate and persist trap results."""
    run = siteUtils.getRunNumber()
    for det_name in det_names:
        file_prefix = '{}_{}'.format(run, det_name)
        trap_file = '%s_traps.fits' % det_name
        eotestUtils.addHeaderData(trap_file, TESTTYPE='TRAP',
                                  DATE=eotestUtils.utc_now_isoformat())

        results.append(siteUtils.make_fileref(trap_file))

        mask_file = '%s_traps_mask.fits' % det_name
        results.append(siteUtils.make_fileref(mask_file))

        results_file = '%s_eotest_results.fits' % det_name
        data = sensorTest.EOTestResults(results_file)
        amps = data['AMP']
        num_traps = data['NUM_TRAPS']

        for amp, ntrap in zip(amps, num_traps):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('traps_BOT'), amp=amp, num_traps=ntrap,
                slot=slot, raft=raft))
    return results


def validate_dark_current(results, det_names):
    """Validate and persist dark current results."""
    run = siteUtils.getRunNumber()
    for det_name in det_names:
        file_prefix = '{}_{}'.format(run, det_name)
        results_file = '%s_eotest_results.fits' % file_prefix
        data = sensorTest.EOTestResults(results_file)

        amps = data['AMP']
        dc95s = data['DARK_CURRENT_95']
        for amp, dc95 in zip(amps, dc95s):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('dark_current_BOT'), amp=amp,
                dark_current_95CL=dc95, slot=slot, raft=raft))

        # Persist the png files.
        metadata = dict(TESTTYPE='DARK', TEST_CATEGORY='EO',
                        DETECTOR=det_name, RUN=run)

        pattern = '{}_*noise*.png'.format(file_prefix)
        results.extend(siteUtils.persist_png_files(pattern, file_prefix,
                                                   metadata=metadata))


if __name__ == '__main__':
    det_names = camera_info.get_det_names()
    results = []
    results = validate_fe55(results, det_names)
    results = validate_read_noise(results, det_names)
    results = validate_bright_defects(results, det_names)
    results = validate_dark_defects(results, det_names)
    results = validate_traps(results, det_names)
    results = validate_dark_current(results, det_names)

    results.extend(siteUtils.jobInfo())

    lcatr.schema.write_file(results)
    lcatr.schema.validate_file()
