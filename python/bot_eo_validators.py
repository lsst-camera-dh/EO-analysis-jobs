"""
Functions to run lcatr.schema validation on BOT EO analysis data products.
"""
from __future__ import print_function
import os
import glob
from collections import OrderedDict
import pickle
import numpy as np
from astropy.io import fits
import lcatr.schema
import siteUtils
import eotestUtils
import lsst.eotest.sensor as sensorTest
from camera_components import camera_info
from .tearing_detection import persist_tearing_png_files
from .bot_eo_analyses import make_file_prefix

__all__ = ['validate_fe55', 'validate_read_noise',
           'validate_bright_defects', 'validate_dark_defects',
           'validate_traps', 'validate_dark_current', 'validate_cte',
           'validate_flat_pairs', 'validate_ptc', 'validate_qe',
           'validate_tearing', 'validate_raft_results']


def report_missing_data(validator, missing_data, components='detectors',
                        total=189):
    """Summarize the missing data for the specified components."""
    if len(missing_data) == total:
        print("{}: missing data for all {} {}".format(validator, total,
                                                      components))
    else:
        print("{}: missing data for {} {}".format(validator, len(missing_data),
                                                  components))
        print(missing_data)


def validate_bias_frames(results, det_names):
    """Validate and persist medianed bias frames."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        file_prefix = make_file_prefix(run, det_name)
        bias_frames = glob.glob('{}_median_bias.fits'.format(file_prefix))
        if bias_frames:
            results.append(siteUtils.make_fileref(bias_frames[0]))
        else:
            missing_det_names.append(det_name)
    report_missing_data('validate_bias_frames', missing_det_names)
    return results


def validate_fe55(results, det_names):
    """Validate and persist fe55 gain and psf results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)

        # The output files from producer script.
        gain_file = '%(file_prefix)s_eotest_results.fits' % locals()
        psf_results_files \
            = glob.glob('%(file_prefix)s_psf_results*.fits' % locals())

        if not os.path.isfile(gain_file) or not psf_results_files:
            # Results for this detector are not available so note
            # that and continue with the others.
            missing_det_names.append(det_name)
            continue
        psf_results = psf_results_files[0]

        rolloff_mask = '%(file_prefix)s_edge_rolloff_mask.fits' % locals()

        output_files = psf_results, rolloff_mask

        # Add/update the metadata to the primary HDU of these files.
        for fitsfile in output_files:
            eotestUtils.addHeaderData(fitsfile, TESTTYPE='FE55',
                                      DATE=eotestUtils.utc_now_isoformat())
        results.extend([lcatr.schema.fileref.make(x) for x in output_files])

        # Persist the png files.
        png_file_list = '{}_fe55_task_png_files.txt'.format(det_name)
        with open(png_file_list, 'r') as input_:
            png_files = [x.strip() for x in input_]
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

    report_missing_data('validate_fe55', missing_det_names)

    return results


def validate_read_noise(results, det_names):
    """Validate and persist read noise results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)

        read_noise_file = '%s_eotest_results.fits' % file_prefix
        if not os.path.isfile(read_noise_file):
            # No data for this detector, so note that and continue
            # with the others.
            missing_det_names.append(det_name)
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
        file_prefix = make_file_prefix(run, raft)
        filename = '%s_overscan_correlations.png' % file_prefix
        results.extend(siteUtils.persist_png_files(filename, file_prefix,
                                                   metadata=metadata))

    report_missing_data("validate_read_noise", missing_det_names)

    return results


def validate_bright_defects(results, det_names):
    """Validate and persist bright defects results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        mask_file = '%s_bright_pixel_mask.fits' % file_prefix
        if not os.path.isfile(mask_file):
            missing_det_names.append(det_name)
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

    report_missing_data("validate_bright_defects", missing_det_names)

    return results

def validate_dark_defects(results, det_names):
    """Validate and persist dark defects results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        mask_file = '%s_dark_pixel_mask.fits' % file_prefix
        if not os.path.isfile(mask_file):
            missing_det_names.append(det_name)
            continue
        eotestUtils.addHeaderData(mask_file, TESTTYPE='SFLAT_500',
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(mask_file))

        superflat = '%s_median_sflat.fits' % file_prefix
        eotestUtils.addHeaderData(superflat,
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(superflat))

        eotest_results = '%s_eotest_results.fits' % file_prefix
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
        metadata = dict(DETECTOR=det_name, RUN=run,
                        TESTTYPE='SFLAT_500', TEST_CATEGORY='EO')
        filename = '%s_superflat_dark_defects.png' % file_prefix
        results.extend(siteUtils.persist_png_files(filename, file_prefix,
                                                   metadata=metadata))

    report_missing_data("validate_dark_defects", missing_det_names)

    return results

def validate_traps(results, det_names):
    """Validate and persist trap results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        trap_file = '%s_traps.fits' % file_prefix
        if not os.path.isfile(trap_file):
            missing_det_names.append(det_name)
            continue
        eotestUtils.addHeaderData(trap_file, TESTTYPE='TRAP',
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(trap_file))

        mask_file = '%s_traps_mask.fits' % file_prefix
        results.append(siteUtils.make_fileref(mask_file))

        results_file = '%s_eotest_results.fits' % file_prefix
        data = sensorTest.EOTestResults(results_file)
        amps = data['AMP']
        num_traps = data['NUM_TRAPS']

        for amp, ntrap in zip(amps, num_traps):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('traps_BOT'), amp=amp, num_traps=ntrap,
                slot=slot, raft=raft))

    report_missing_data("validate_traps", missing_det_names)

    return results


def validate_dark_current(results, det_names):
    """Validate and persist dark current results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        results_file = '%s_eotest_results.fits' % file_prefix
        if not os.path.isfile(results_file):
            missing_det_names.append(det_name)
            continue
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

        pattern = '{}_noise.png'.format(file_prefix)
        results.extend(siteUtils.persist_png_files(pattern, file_prefix,
                                                   metadata=metadata))
        pattern = '{}_total_noise_hists.png'.format(file_prefix)
        results.extend(siteUtils.persist_png_files(pattern, file_prefix,
                                                   metadata=metadata))

    report_missing_data("validate_dark_current", missing_det_names)

    return results


def validate_cte(results, det_names):
    """Validate the CTE task results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        superflats \
            = sorted(glob.glob('{}_superflat_*.fits'.format(file_prefix)))
        if not superflats:
            missing_det_names.append(det_name)
            continue
        for item in superflats:
            eotestUtils.addHeaderData(item, FILENAME=item,
                                      DATE=eotestUtils.utc_now_isoformat())
        results.extend([siteUtils.make_fileref(x) for x in superflats])

        results_file = '%s_eotest_results.fits' % file_prefix
        data = sensorTest.EOTestResults(results_file)
        amps = data['AMP']

        cti_high_serial = data['CTI_HIGH_SERIAL']
        cti_high_serial_error = data['CTI_HIGH_SERIAL_ERROR']
        cti_high_parallel = data['CTI_HIGH_PARALLEL']
        cti_high_parallel_error = data['CTI_HIGH_PARALLEL_ERROR']

        cti_low_serial = data['CTI_LOW_SERIAL']
        cti_low_serial_error = data['CTI_LOW_SERIAL_ERROR']
        cti_low_parallel = data['CTI_LOW_PARALLEL']
        cti_low_parallel_error = data['CTI_LOW_PARALLEL_ERROR']

        for values in zip(amps,
                          cti_high_serial, cti_high_serial_error,
                          cti_high_parallel, cti_high_parallel_error,
                          cti_low_serial, cti_low_serial_error,
                          cti_low_parallel, cti_low_parallel_error):
            results.append(lcatr.schema.valid(lcatr.schema.get('cte_BOT'),
                                              amp=values[0],
                                              cti_high_serial=values[1],
                                              cti_high_serial_error=values[2],
                                              cti_high_parallel=values[3],
                                              cti_high_parallel_error=values[4],
                                              cti_low_serial=values[5],
                                              cti_low_serial_error=values[6],
                                              cti_low_parallel=values[7],
                                              cti_low_parallel_error=values[8],
                                              slot=slot, raft=raft))

        # Persist the png files.
        png_file_list = '{}_cte_task_png_files.txt'.format(det_name)
        with open(png_file_list, 'r') as input_:
            png_files = [x.strip() for x in input_]
        metadata = dict(DETECTOR=det_name, RUN=run,
                        TESTTYPE='SFLAT_500', TEST_CATEGORY='EO')
        results.extend(siteUtils.persist_png_files('', file_prefix,
                                                   png_files=png_files,
                                                   metadata=metadata))

    report_missing_data("validate_cte", missing_det_names)

    return results


def validate_flat_pairs(results, det_names):
    """Validate the flat pair analysis results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        det_resp_data = '%s_det_response.fits' % file_prefix
        if not os.path.isfile(det_resp_data):
            missing_det_names.append(det_name)
            continue
        eotestUtils.addHeaderData(det_resp_data, DETECTOR=det_name,
                                  TESTTYPE='FLAT',
                                  DATE=eotestUtils.utc_now_isoformat())
        results.append(siteUtils.make_fileref(det_resp_data))

        results_file = '%s_eotest_results.fits' % file_prefix
        data = sensorTest.EOTestResults(results_file)
        amps = data['AMP']
        full_well_data = data['FULL_WELL']
        max_frac_dev_data = data['MAX_FRAC_DEV']

        for amp, full_well, max_frac_dev in zip(amps, full_well_data,
                                                max_frac_dev_data):
            results.append(lcatr.schema.valid(
                lcatr.schema.get('flat_pairs_BOT'),
                amp=amp, full_well=full_well, max_frac_dev=max_frac_dev,
                slot=slot, raft=raft))

        # Persist the png files.
        metadata = dict(DETECTOR=det_name, RUN=run,
                        TESTTYPE='FLAT', TEST_CATEGORY='EO')
        results.extend(siteUtils.persist_png_files(('%s_linearity*.png'
                                                    % file_prefix),
                                                   file_prefix,
                                                   metadata=metadata))

    report_missing_data("validate_flat_pairs", missing_det_names)

    return results


def validate_ptc(results, det_names):
    """Validate the PTC results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)
        ptc_results = '%s_ptc.fits' % file_prefix
        if not os.path.isfile(ptc_results):
            missing_det_names.append(det_name)
            continue
        eotestUtils.addHeaderData(ptc_results, TESTTYPE='FLAT',
                                  DATE=eotestUtils.utc_now_isoformat())

        results.append(siteUtils.make_fileref(ptc_results))

        results_file = '%s_eotest_results.fits' % file_prefix
        data = sensorTest.EOTestResults(results_file)

        columns = (data['AMP'], data['PTC_GAIN'], data['PTC_GAIN_ERROR'],
                   data['PTC_A00'], data['PTC_A00_ERROR'], data['PTC_NOISE'],
                   data['PTC_NOISE_ERROR'], data['PTC_TURNOFF'])
        for amp, gain, gain_error, a00, a00_error,\
            noise, noise_error, turnoff in zip(*columns):
            results.append(lcatr.schema.valid(lcatr.schema.get('ptc_BOT'),
                                              amp=amp, ptc_gain=gain,
                                              ptc_gain_error=gain_error,
                                              ptc_a00=a00,
                                              ptc_a00_error=a00_error,
                                              ptc_noise=noise,
                                              ptc_noise_error=noise_error,
                                              ptc_turnoff=turnoff,
                                              slot=slot, raft=raft))
        # Persist the png files.
        metadata = dict(DETECTOR=det_name, RUN=run,
                        TESTTYPE='FLAT', TEST_CATEGORY='EO')

        results.extend(siteUtils.persist_png_files('%s*ptcs.png' % file_prefix,
                                                   file_prefix,
                                                   metadata=metadata))

    report_missing_data("validate_ptc", missing_det_names)

    return results


def validate_qe(results, det_names):
    """Validate the QE results."""
    run = siteUtils.getRunNumber()
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)

        qe_results_file = '%s_QE.fits' % file_prefix
        if not os.path.isfile(qe_results_file):
            missing_det_names.append(det_name)
            continue
        with fits.open(qe_results_file) as qe_results:
            qe_data = qe_results['QE_BANDS'].data
            QE = OrderedDict((band, []) for band in qe_data.field('BAND'))
            for amp in range(1, 17):
                values = qe_data.field('AMP%02i' % amp)
            for band, value in zip(QE, values):
                QE[band].append(value)

        for band in QE:
            for amp in range(1, 17):
                results.append(lcatr.schema.valid(
                    lcatr.schema.get('qe_BOT_analysis'),
                    band=band, QE=QE[band][amp-1],
                    amp=amp, slot=slot, raft=raft))

        qe_files = glob.glob('%s_*QE*.fits' % file_prefix)
        for item in qe_files:
            eotestUtils.addHeaderData(item, TESTTYPE='LAMBDA',
                                      DATE=eotestUtils.utc_now_isoformat())
        results.extend([siteUtils.make_fileref(item) for item in qe_files])

        # Persist the png files.
        metadata = dict(DETECTOR=det_name, RUN=run,
                        TESTTYPE='LAMBDA', TEST_CATEGORY='EO')
        results.extend(siteUtils.persist_png_files('%s*qe.png' % file_prefix,
                                                   file_prefix,
                                                   metadata=metadata))
        results.extend(siteUtils.persist_png_files('%s*flat.png' % file_prefix,
                                                   file_prefix,
                                                   metadata=metadata))

    report_missing_data("validate_qe", missing_det_names)

    return results


def validate_tearing(results, det_names):
    """Validate the tearing analysis results."""
    run = siteUtils.getRunNumber()
    schema = lcatr.schema.get('tearing_detection_BOT')
    missing_det_names = []
    for det_name in det_names:
        raft, slot = det_name.split('_')
        file_prefix = make_file_prefix(run, det_name)

        tearing_results_file = '%s_tearing_stats.pkl' % file_prefix
        if not os.path.isfile(tearing_results_file):
            missing_det_names.append(det_name)
            continue
        with open(tearing_results_file, 'rb') as input_:
            tearing_stats = pickle.load(input_)
        for values in tearing_stats:
            stats = dict(kv for kv in zip(('job_name', 'subset', 'sensor_id',
                                           'detections', 'slot', 'raft'),
                                          list(values) + [slot, raft]))
            results.append(lcatr.schema.valid(schema, **stats))

    png_files = sorted(glob.glob('*_tearing.png'))
    results.extend(persist_tearing_png_files(png_files))

    report_missing_data("validate_tearing", missing_det_names)

    return results


def validate_raft_results(results, raft_names):
    """Validate the raft level results."""
    run = siteUtils.getRunNumber()
    slot_names = camera_info.get_slot_names()
    md = siteUtils.DataCatalogMetadata(ORIGIN=siteUtils.getSiteName(),
                                       TEST_CATEGORY='EO',
                                       DATA_PRODUCT='EOTEST_RESULTS')
    missing_raft_names = []
    for raft_name in raft_names:
        for slot_name in slot_names:
            det_name = make_file_prefix(raft_name, slot_name)
            file_prefix = make_file_prefix(run, det_name)
            results_file = '{}_eotest_results.fits'.format(file_prefix)
            if not os.path.isfile(results_file):
                if raft_name not in missing_raft_names:
                    missing_raft_names.append(raft_name)
                continue
            eotestUtils.addHeaderData(results_file, DETECTOR=det_name,
                                      DATE=eotestUtils.utc_now_isoformat(),
                                      RUNNUM=run)
            results.append(lcatr.schema.fileref.make(
                results_file, metadata=md(SLOT=slot_name, RAFT=raft_name)))

        # Persist the png files.
        png_file_list = '{}_raft_results_task_png_files.txt'.format(raft_name)
        if not os.path.isfile(png_file_list):
            continue
        with open(png_file_list, 'r') as input_:
            png_files = [x.strip() for x in input_]
        metadata = dict(TEST_CATEGORY='EO', DETECTOR=det_name, RUN=run)
        results.extend(siteUtils.persist_png_files('', file_prefix,
                                                   png_files=png_files,
                                                   metadata=metadata))

    report_missing_data("validate_raft_results", missing_raft_names,
                        components='rafts', total=21)

    return results


if __name__ == '__main__':
    det_names = camera_info.get_det_names()
    raft_names = camera_info.get_raft_names()

    results = []
    results = validate_fe55(results, det_names)
    results = validate_read_noise(results, det_names)
    results = validate_bright_defects(results, det_names)
    results = validate_dark_defects(results, det_names)
    results = validate_traps(results, det_names)
    results = validate_dark_current(results, det_names)
    results = validate_cte(results, det_names)
    results = validate_flat_pairs(results, det_names)
    results = validate_ptc(results, det_names)
    results = validate_qe(results, det_names)
    results = validate_tearing(results, det_names)
    results = validate_raft_results(results, raft_names)

    results.extend(siteUtils.jobInfo())

    lcatr.schema.write_file(results)
    lcatr.schema.validate_file()
