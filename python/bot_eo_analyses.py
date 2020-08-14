"""
Producer script for BOT analyses.
"""
import os
import re
import glob
import copy
import time
import shutil
import json
import pickle
import warnings
from collections import defaultdict
import configparser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
from lsst.eotest.sensor.EOTestPlots import OverscanTestPlots
import eotestUtils
import siteUtils
from correlated_noise import correlated_noise, raft_level_oscan_correlations
from camera_components import camera_info
from tearing_detection import tearing_detection
from multiprocessor_execution import run_device_analysis_pool
try:
    import scope
    import multiscope
except ImportError:
    print("scope and/or multiscope not imported")

__all__ = ['make_file_prefix',
           'glob_pattern',
           'get_mask_files',
           'get_amplifier_gains',
           'medianed_dark_frame',
           'bias_filename',
           'get_raft_files_by_slot',
           'get_analysis_types',
           'fe55_task',
           'bias_frame_task',
           'bias_stability_task',
           'scan_mode_analysis_task',
           'get_scan_mode_files',
           'read_noise_task',
           'raft_noise_correlations',
           'bright_defects_task',
           'dark_defects_task',
           'traps_task',
           'dark_current_task',
           'plot_ccd_total_noise',
           'cte_task',
           'plot_cte_results',
           'find_flat2_bot',
           'flat_pairs_task',
           'nonlinearity_task',
           'ptc_task',
           'bf_task',
           'qe_task',
           'tearing_task',
           'overscan_task',
           'persistence_task',
           'repackage_summary_files',
           'mondiode_value',
           'get_nlc_func',
           'run_jh_tasks',
           'run_python_task_or_cl_script']


class GlobPattern:
    """Functor class to return glob patterns for finding BOT_acq data."""
    def __init__(self):
        config = configparser.ConfigParser(inline_comment_prefixes=("#", ))
        config.optionxform = str
        default_config = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                                      'data', 'BOT_jh_glob_patterns.ini')
        cfg_file = os.environ.get('LCATR_JH_GLOB_PATTERN_FILE', default_config)
        config.read(cfg_file)
        self.task_patterns = dict(config.items('BOT_acq'))

    def __call__(self, task, det_name):
        pattern = self.task_patterns[task]
        return '{}/*_{}.fits'.format(pattern, det_name)


glob_pattern = GlobPattern()


def mondiode_value(flat_file, exptime, factor=5,
                   pd_filename='Photodiode_Readings*.txt'):
    """
    Compute the mean current measured by the monitoring photodiode.

    Parameters
    ----------
    flat_file: str
        Path to the flat frame FITS file.   The pd data file
        is assumed to be in the same directory.
    exptime: float
        Exposure time in seconds.
    factor: float [5]
        Factor to use to extract the baseline current values from the
        data using ythresh = (min(y) + max(y))/factor + min(y)
    pd_filename: str ['Photodiode_Readings*.txt']
        Basename of photodiode readings file.

    Returns
    -------
    float: The mean current.
    """
    # Try reading the MONDIODE keyword first.
    with fits.open(flat_file) as hdulist:
        if 'MONDIODE' in hdulist[0].header:
            return hdulist[0].header['MONDIODE']

    # Compute the value from the photodiode readings file.
    pd_file = glob.glob(os.path.join(os.path.dirname(flat_file),
                                     pd_filename))[0]
    x, y = np.recfromtxt(pd_file).transpose()
    # Threshold for finding baseline current values:
    ythresh = (min(y) + max(y))/factor + min(y)
    # Subtract the median of the baseline values to get a calibrated
    # current.
    y -= np.median(y[np.where(y < ythresh)])
    integral = sum((y[1:] + y[:-1])/2*(x[1:] - x[:-1]))
    return integral/exptime


def get_mask_files(det_name):
    """
    Get the mask files from previous jobs for the specified sensor.

    Parameters
    ----------
    det_name: str
        The detector name in the focal plane, e.g., 'R22_S11'.

    Returns
    -------
    list of mask file paths.
    """
    badpixel_run = siteUtils.get_analysis_run('badpixel')
    bias_run = siteUtils.get_analysis_run('bias')

    if badpixel_run is not None or bias_run is not None:
        with open('hj_fp_server.pkl', 'rb') as fd:
            hj_fp_server = pickle.load(fd)

    if badpixel_run is not None:
        mask_files = hj_fp_server.get_files('pixel_defects_BOT',
                                            f'{det_name}*mask*.fits',
                                            run=badpixel_run)
        print(f"Mask files from run {badpixel_run} and {det_name}:")
        for item in mask_files:
            print(item)
        print()
    else:
        mask_files = siteUtils.dependency_glob(f'*{det_name}*mask*.fits',
                                               jobname='pixel_defects_BOT',
                                               description='pixel defects masks')
    if bias_run is not None:
        rolloff_mask_files = hj_fp_server.get_files('bias_frame_BOT',
                                                    f'{det_name}_*mask*.fits',
                                                    run=bias_run)
        print(f"Edge rolloff mask file from run {bias_run} and {det_name}:")
        for item in rolloff_mask_files:
            print(item)
        print()
    else:
        rolloff_mask_files = siteUtils.dependency_glob(f'*{det_name}*mask*.fits',
                                                       description='rolloff masks',
                                                       jobname='bias_frame_BOT')
    mask_files.extend(rolloff_mask_files)

    return mask_files


def make_file_prefix(run, component_name):
    """
    Compose the run number and component name into string prefix
    to use with filenames.
    """
    return "{}_{}".format(component_name, run)


def make_bias_filename(run, det_name):
    """Make the bias filename from the run and det_name."""
    file_prefix = make_file_prefix(run, det_name)
    return f'{file_prefix}_median_bias.fits'


def medianed_dark_frame(det_name):
    """
    The medianed dark frame from the pixel defects task.
    """
    pattern = f'{det_name}*_median_dark_bp.fits'
    dark_run = siteUtils.get_analysis_run('dark')
    if dark_run is None:
        return siteUtils.dependency_glob(pattern, description='Dark frame:')[0]

    # Retrieve bias file from previous run.
    with open('hj_fp_server.pkl', 'rb') as fd:
        hj_fp_server = pickle.load(fd)
    filename = hj_fp_server.get_files('pixel_defects_BOT', pattern,
                                      run=dark_run)[0]
    print("Dark frame:")
    print(filename)
    return filename


def bias_filename(run, det_name):
    """
    The bias frame file derived from stacked bias files.
    """
    bias_run = siteUtils.get_analysis_run('bias')
    if bias_run is None:
        filename = make_bias_filename(run, det_name)
        if not os.path.isfile(filename):
            # Look for bias file from prerequisite job.
            return siteUtils.dependency_glob(filename,
                                             description='Bias frames:')[0]
    else:
        # Retrieve bias file from previous run.
        with open('hj_fp_server.pkl', 'rb') as fd:
            hj_fp_server = pickle.load(fd)
        filename = hj_fp_server.get_files('bias_frame_BOT',
                                          f'*{det_name}*median_bias.fits',
                                          run=bias_run)[0]
    print("Bias frame:")
    print(filename)
    return filename


def fe55_task(run, det_name, fe55_files):
    "Single sensor execution of the Fe55 analysis task."
    file_prefix = make_file_prefix(run, det_name)
    title = '{}, {}'.format(run, det_name)

    bias_frame = bias_filename(run, det_name)
    png_files = []

    try:
        pixel_stats = sensorTest.Fe55PixelStats(fe55_files,
                                                sensor_id=file_prefix)
        png_files.append('%s_fe55_p3_p5_hists.png' % file_prefix)
        siteUtils.make_png_file(pixel_stats.pixel_hists, png_files[-1],
                                pix0='p3', pix1='p5')

        png_files.append('%s_fe55_p3_p5_profiles.png' % file_prefix)
        siteUtils.make_png_file(pixel_stats.pixel_diff_profile,
                                png_files[-1], pixel_coord='x',
                                pix0='p3', pix1='p5')

    except Exception:
        # Encountered error processing data or generating pngs so skip
        # these plots.
        pass

    rolloff_mask_file = '%s_edge_rolloff_mask.fits' % file_prefix
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    hist_nsig = 20

    task = sensorTest.Fe55Task()
    task.config.temp_set_point = -100.
    task.run(file_prefix, fe55_files, (rolloff_mask_file,),
             bias_frame=bias_frame, accuracy_req=0.01, hist_nsig=hist_nsig,
             linearity_correction=get_nlc_func(det_name))

    # Fe55 gain and psf analysis results plots for the test report.
    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)

    png_files.append('%s_gains.png' % file_prefix)
    siteUtils.make_png_file(plots.gains, png_files[-1])

    png_files.append('%s_fe55_median_bias.png' % file_prefix)
    siteUtils.make_png_file(sensorTest.plot_flat, png_files[-1], bias_frame,
                            title='%s, median bias frame' % title,
                            annotation='ADU/pixel, overscan-subtracted')

    fe55_file = glob.glob('%s_psf_results*.fits' % file_prefix)[0]
    png_files.append('%s_fe55_dists.png' % file_prefix)
    siteUtils.make_png_file(plots.fe55_dists, png_files[-1],
                            fe55_file=fe55_file, xrange_scale=5,
                            hist_nsig=hist_nsig)

    png_files.append('%s_psf_dists.png' % file_prefix)
    siteUtils.make_png_file(plots.psf_dists, png_files[-1],
                            fe55_file=fe55_file)

    png_file_list = '{}_fe55_task_png_files.txt'.format(det_name)
    with open(png_file_list, 'w') as output:
        for item in png_files:
            if os.path.isfile(item):
                output.write('{}\n'.format(item))


def flat_gain_stability_task(run, det_name, flat_files, mask_files=(),
                             bias_frame=None, dark_frame=None,
                             mondiode_func=mondiode_value):
    """
    Task to compute the median signal level in each amp of each flat
    frame as a function of mjd and seqnum.

    Parameters
    ----------
    run: str
        Run number.
    det_name: str
        Sensor name in the focal plane, e.g., 'R22_S11'.
    flat_files: list
        Sequence of flat files all taken with the same exposure time and
        incident flux.
    mask_files: list-like [()]
        Mask files to apply for computing image statistics.
    bias_frame: str [None]
        Medianed bias frame to use for bias subtraction.
    dark_frame: str [None]
        Medianed dark frame to use for dark subtraction.
    mondiode_func: function [bot_eo_analyses.mondiode_value]
        Function to use for computing monitoring diode current.

    Returns
    -------
    (pandas.DataFrame, output_file) containing the time history of the median
        signals from each amp and the name of the output pickle file containing
        these data.
    """
    file_prefix = make_file_prefix(run, det_name)
    outfile = f'{file_prefix}_flat_signal_sequence.pickle'
    df = sensorTest.flat_signal_sequence(flat_files, bias_frame=bias_frame,
                                         dark_frame=dark_frame,
                                         mask_files=mask_files,
                                         mondiode_func=mondiode_func)
    df.to_pickle(outfile)
    return df, outfile


def gain_stability_task(run, det_name, fe55_files):
    """
    This task fits the Fe55 clusters to the cluster data from each frame
    sequence and writes a pickle file with the gains as a function of
    sequence number and MJD-OBS.

    Parameters
    ----------
    run: str
        Run number.
    det_name: str
        Sensor name in the focal plane, e.g., 'R22_S11'.
    fe55_files: list
        Raw Fe55 for the sensor being consider.  The MJD-OBS values
        will be extracted from these files.

    Returns:
    (pandas.DataFrame, str), i.e., a tuple of the data frame containing
    the gain sequence and the file name of the output pickle file.
    """
    file_prefix = make_file_prefix(run, det_name)

    # Extract MJD-OBS values into a dict to provide look up table in
    # case there are missing sequence frames in the psf results table.
    mjd_obs = dict()
    for item in fe55_files:
        with fits.open(item) as hdus:
            mjd_obs[hdus[0].header['SEQNUM']] = hdus[0].header['MJD-OBS']

    psf_results_file = sorted(glob.glob(f'{file_prefix}_psf_results*.fits'))[0]
    df = sensorTest.gain_sequence(det_name, psf_results_file)
    df['mjd'] = [mjd_obs[seqnum] for seqnum in df['seqnum']]
    outfile = f'{file_prefix}_gain_sequence.pickle'
    df.to_pickle(outfile)

    return df, outfile


class GetAmplifierGains:
    """
    Functor class to provide gains either from an ET db lookup or
    from the EO test results file, depending on the traveler instance
    configuration.
    """
    def __init__(self, bot_eo_config_file=None,
                 et_results_file='et_results.pkl', run=None):
        if run is None:
            self.run = siteUtils\
                .get_analysis_run('gain', bot_eo_config_file=bot_eo_config_file)
        else:
            self.run = run

        if self.run is not None:
            if self.run.endswith('.json'):
                self._get_curated_gains()
            else:
                self.curated_gains = None
                self._get_gains_from_run(et_results_file)

    def _get_curated_gains(self):
        gain_file = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                                 'data', self.run)
        # Read the curated gains from the json file.
        with open(gain_file, 'r') as fd:
            self.curated_gains = json.load(fd)
        print("GetAmplifierGains: Using gains from", gain_file)

    def _get_gains_from_run(self, et_results_file):
        if not os.path.isfile(et_results_file):
            self.et_results = siteUtils.ETResults(self.run)
            with open(et_results_file, 'wb') as fd:
                pickle.dump(self.et_results, fd)
        else:
            with open(et_results_file, 'rb') as fd:
                self.et_results = pickle.load(fd)
        print("GetAmplifierGains: Using gains from run", self.run)

    def __call__(self, file_pattern):
        if self.run is None:
            return _get_amplifier_gains(file_pattern)
        # Extract the det_name from the file pattern.
        match = re.search('R\d\d_S\w\d', file_pattern)
        if match is None:
            message = f"no det_name match in {file_pattern}"
            raise RuntimeError("GetAmplifierGains.__call__: " + message)
        det_name = file_pattern[match.start(): match.end()]

        if self.curated_gains is not None:
            print('GetAmplifierGains.__call__: retrieving curated gains.')
            my_gains = self.curated_gains[det_name]
            if len(my_gains) == 8:
                channels = siteUtils.ETResults.wf_amp_names
            else:
                channels = siteUtils.ETResults.amp_names
            return {amp: my_gains[_] for amp, _ in enumerate(channels, 1)}

        print("GetAmplifierGains.__call__: retrieving gains from eT.")
        gains = self.et_results.get_amp_gains(det_name)
        if not gains:
            print("GetAmplifierGains.__call__: "
                  "gains from eT not found, using unit gains.")
            return {amp: 1 for amp in range(1, 17)}
        return gains


def _get_amplifier_gains(file_pattern=None):
    """Extract the gains for each amp in an eotest_results file."""
    if (os.environ.get('LCATR_USE_UNIT_GAINS', 'False') == 'True'
        or file_pattern is None):
        print("_get_amplifier_gains: using unit gains")
        return {amp: 1 for amp in range(1, 17)}

    # Attempt to retrieve gains from fe55_analysis_BOT then ptc_BOT.
    # If neither are available, then use unit gains.
    print("_get_amplifier_gains: trying fe55_analysis_BOT")
    results_files = siteUtils.dependency_glob(file_pattern,
                                              jobname='fe55_analysis_BOT')
    if not results_files:
        print("_get_amplifier_gains: trying ptc_BOT")
        results_files = siteUtils.dependency_glob(file_pattern,
                                                  jobname='ptc_BOT')
    if not results_files:
        print("_get_amplifier_gains: both fe55 and ptc retrievals failed. "
              "using unit gains.")
        return {amp: 1 for amp in range(1, 17)}

    eotest_results_file = results_files[0]
    data = sensorTest.EOTestResults(eotest_results_file)
    amps = data['AMP']
    gains = data['GAIN']
    return dict(zip(amps, gains))

try:
    get_amplifier_gains = GetAmplifierGains()
except KeyError as eobj:
    warnings.warn(f"KeyError: {eobj} when creating get_amplifier_gains "
                  "function object.  Assigning to default function.")
    get_amplifier_gains = _get_amplifier_gains


def bias_frame_task(run, det_name, bias_files, bias_frame=None):
    """Create a median bias file for use by downstream tasks."""
    if bias_frame is None:
        bias_frame = make_bias_filename(run, det_name)
    amp_geom = sensorTest.makeAmplifierGeometry(bias_files[0])
    imutils.superbias_file(bias_files, amp_geom.serial_overscan, bias_frame)
    file_prefix = make_file_prefix(run, det_name)
    rolloff_mask_file = f'{file_prefix}_edge_rolloff_mask.fits'
    sensorTest.rolloff_mask(bias_files[0], rolloff_mask_file)

def image_stats(image, nsigma=10):
    """Compute clipped mean and stdev of the image."""
    stat_ctrl = afwMath.StatisticsControl(numSigmaClip=nsigma)
    flags = afwMath.MEANCLIP | afwMath.STDEVCLIP
    stats = afwMath.makeStatistics(image, flags=flags, sctrl=stat_ctrl)
    return stats.getValue(afwMath.MEANCLIP), stats.getValue(afwMath.STDEVCLIP)

def bias_stability_task(run, det_name, bias_files, nsigma=10):
    """Compute amp-wise bias stability time histories and serial profiles."""
    raft, slot = det_name.split('_')
    file_prefix = make_file_prefix(run, det_name)
    data = defaultdict(list)

    fig = plt.figure(figsize=(16, 16))
    xlabel_amps = (13, 14, 15, 16)
    ylabel_amps = (1, 5, 9, 13)
    ax = {amp: fig.add_subplot(4, 4, amp) for amp in range(1, 17)}

    for bias_file in bias_files:
        ccd = sensorTest.MaskedCCD(bias_file)
        for amp in ccd:
            # Retrieve the per row overscan subtracted imaging section.
            amp_image = ccd.unbiased_and_trimmed_image(amp)
            # Plot the median of each column versus serial pixel number.
            imarr = amp_image.getImage().array
            ax[amp].plot(range(imarr.shape[1]), np.median(imarr, axis=0))
            # Compute 10-sigma clipped mean and stdev
            mean, stdev = image_stats(amp_image, nsigma=nsigma)
            data['raft'].append(raft)
            data['slot'].append(slot)
            data['tseqnum'].append(ccd.md.get('TSEQNUM'))
            data['amp'].append(amp)
            data['mean'].append(mean)
            data['stdev'].append(stdev)
    plt.suptitle(f'{det_name}, Run {run}\nmedian signal (ADU) vs column')
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    for amp in ccd:
        ax[amp].annotate(f'amp {amp}', (0.5, 0.95),
                         xycoords='axes fraction', ha='center')
    plt.savefig(f'{file_prefix}_bias_serial_profiles.png')
    df = pd.DataFrame(data=data)
    df.to_pickle(f'{file_prefix}_bias_frame_stats.pkl')


def get_scan_mode_files(raft_name):
    """Get the scan mode filenames for a given raft, organized by
    scan mode folder and slot.
    """
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    pattern = glob_pattern('scan_mode', '{}_*'.format(raft_name))
    files = siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname)
    scan_mode_files = defaultdict(dict)
    for item in files:
        scan_dir = os.path.basename(os.path.dirname(item))
        slot = os.path.basename(item).split('.')[0].split('_')[-1]
        scan_mode_files[scan_dir][slot] = item
    return scan_mode_files


def get_raft_arrays(raft_files):
    """Get raftarrrays list for passing to scan mode plotting code."""
    seglist = multiscope.slot_ids()
    raftarrays = []
    for slot in ['S{}'.format(seg) for seg in seglist]:
        raftarrays.append(scope.get_scandata_fromfile(raft_files[slot]))
    return raftarrays, seglist


def scan_mode_analysis_task(run, raft_name, scan_mode_files):

    """Scan mode analysis task."""
    file_prefix = make_file_prefix(run, raft_name)
    for scan_dir, raft_files in scan_mode_files.items():
        raft_arrays, seg_list = get_raft_arrays(raft_files)
        # Make the dispersion plots, one per sensor.
        for seg, scandata in zip(seg_list, raft_arrays):
            slot = 'S' + seg
            det_name = '_'.join((raft_name, slot))
            disp_plot_title = '{}, Run {}, {}'.format(det_name, run, scan_dir)
            scope.plot_scan_dispersion(scandata, title=disp_plot_title)
            disp_outfile \
                = '{}_{}_{}_dispersion.png'.format(det_name, run, scan_dir)
            plt.savefig(disp_outfile)
            plt.close()
        # Make the multiscope plots for each raft.
        title = '{}, Run {}, {}'.format(raft_name, run, scan_dir)
        multiscope.plot_raft_allchans(raft_arrays, seg_list, suptitle=title)
        outfile = '{}_{}_multiscope.png'.format(file_prefix, scan_dir)
        plt.savefig(outfile)
        plt.close()


def read_noise_task(run, det_name, bias_files, gains, mask_files=(),
                    system_noise=None):
    """Run the read noise tasks on a single detector."""
    file_prefix = make_file_prefix(run, det_name)
    title = '{}, {}'.format(run, det_name)

    task = sensorTest.ReadNoiseTask()
    task.config.temp_set_point = -100.
    task.run(file_prefix, bias_files, gains, system_noise=system_noise,
             mask_files=mask_files, use_overscan=True)

    # Compute amp-amp correlated noise.
    my_bias_files = copy.deepcopy(bias_files)
    _, corr_fig, _ = correlated_noise(my_bias_files, target=0,
                                      make_plots=True, title=title)
    plt.figure(corr_fig.number)
    plt.savefig('%s_correlated_noise.png' % file_prefix)


def raft_noise_correlations(run, raft_name, bias_file_dict):
    """Raft-level noise-correlation analysis."""
    file_prefix = make_file_prefix(run, raft_name)
    title = "Overscan correlations, Run {}, {}".format(run, raft_name)
    raft_level_oscan_correlations(bias_file_dict, title=title)
    plt.savefig('{}_overscan_correlations.png'.format(file_prefix))


def bright_defects_task(run, det_name, dark_files, gains, mask_files=(),
                        bias_frame=None):
    "Single sensor execution of the bright pixels task."
    file_prefix = make_file_prefix(run, det_name)
    title = '{}, {}'.format(run, det_name)

    task = sensorTest.BrightPixelsTask()
    task.config.temp_set_point = -100.
    task.run(file_prefix, dark_files, mask_files, gains, bias_frame=bias_frame,
             linearity_correction=get_nlc_func(det_name))

    title = '%s, medianed dark for bright defects analysis' % file_prefix
    annotation = 'e-/pixel, gain-corrected, bias-subtracted'
    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_medianed_dark.png' % file_prefix,
                            '%s_median_dark_bp.fits' % file_prefix,
                            title=title, annotation=annotation,
                            bias_frame=bias_frame, gains=gains, binsize=4)


def dark_defects_task(run, det_name, sflat_files, mask_files=(),
                      bias_frame=None):
    """Single sensor execution of the dark defects task."""
    file_prefix = make_file_prefix(run, det_name)
    title = '{}, {}'.format(run, det_name)

    task = sensorTest.DarkPixelsTask()
    task.run(file_prefix, sflat_files, mask_files, bias_frame=bias_frame,
             linearity_correction=get_nlc_func(det_name))

    title = '%s, superflat for dark defects analysis' % file_prefix
    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_superflat_dark_defects.png' % file_prefix,
                            '%s_median_sflat.fits' % file_prefix,
                            title=title, annotation='ADU/pixel',
                            flatten=True, binsize=4)

def traps_task(run, det_name, trap_file, gains, mask_files=(), bias_frame=None):
    """Single sensor execution of the traps analysis task."""
    file_prefix = make_file_prefix(run, det_name)
    task = sensorTest.TrapTask()
    task.run(file_prefix, trap_file, mask_files, gains, bias_frame=bias_frame,
             linearity_correction=get_nlc_func(det_name))


def dark_current_task(run, det_name, dark_files, gains, mask_files=(),
                      temp_set_point=-100., bias_frame=None):
    """Single sensor execution of the dark current task."""
    file_prefix = make_file_prefix(run, det_name)
    task = sensorTest.DarkCurrentTask()
    task.config.temp_set_point = temp_set_point
    return task.run(file_prefix, dark_files, mask_files, gains,
                    bias_frame=bias_frame,
                    linearity_correction=get_nlc_func(det_name))


def plot_ccd_total_noise(run, det_name, dark_curr_pixels, dark95s,
                         eotest_results_file):
    """
    Make CCD-level total noise summary plots using the dark current
    measurements and an existing eotest results file containing
    the read noise measurements.
    """
    file_prefix = make_file_prefix(run, det_name)
    plots = sensorTest.EOTestPlots(det_name, results_file=eotest_results_file)
    siteUtils.make_png_file(plots.total_noise, '%s_noise.png' % file_prefix,
                            dark95s=dark95s)


def cte_task(run, det_name, sflat_files, gains, mask_files=(),
             flux_level='high', bias_frame=None):
    """Single sensor execution of the CTE task."""
    file_prefix = make_file_prefix(run, det_name)

    task = sensorTest.CteTask()
    task.run(file_prefix, sflat_files, flux_level=flux_level, gains=gains,
             mask_files=mask_files, bias_frame=bias_frame,
             linearity_correction=get_nlc_func(det_name))
    # TODO: refactor CteTask to provide access to the superflat filename
    # instead of recomputing it here.
    superflat_file = '{}_superflat_{}.fits'.format(file_prefix, flux_level)
    return superflat_file


def plot_cte_results(run, det_name, superflat_file, eotest_results_file,
                     mask_files=()):
    """
    Create single CCD superflat mosaic plots and plots of the serial and
    parallel CTE profiles.
    """
    file_prefix = make_file_prefix(run, det_name)
    flux_level = 'low' if 'low' in os.path.basename(superflat_file) else 'high'
    plots \
        = sensorTest.EOTestPlots(file_prefix, results_file=eotest_results_file)

    png_files = []
    png_files.append(superflat_file.replace('.fits', '.png'))
    siteUtils.make_png_file(sensorTest.plot_flat, png_files[-1],
                            superflat_file,
                            title=('%s, %s, CTE superflat, %s flux '
                                   % (run, det_name, flux_level)),
                            annotation='ADU/pixel', flatten=True, binsize=4)

    png_files.append('%s_serial_oscan_%s.png' % (file_prefix, flux_level))
    siteUtils.make_png_file(plots.cte_profiles, png_files[-1], flux_level,
                            superflat_file, mask_files, serial=True)

    png_files.append('%s_parallel_oscan_%s.png' % (file_prefix, flux_level))
    siteUtils.make_png_file(plots.cte_profiles, png_files[-1], flux_level,
                            superflat_file, mask_files, serial=False)

    return png_files


def find_flat2_bot(file1):
    """
    Function to find the flat2 file corresponding to the flat1 frame
    for BOT-level acquisitions.
    """
    basename_pattern = '*_' + file1[-len('R22_S11.fits'):]
    pattern = os.path.join(file1.split('flat1')[0] + 'flat0*',
                           basename_pattern)
    try:
        flat2 = glob.glob(pattern)[0]
    except IndexError:
        pattern = os.path.join(file1.split('flat1')[0] + 'flat2*',
                               basename_pattern)
        flat2 = glob.glob(pattern)[0]
    return flat2


def flat_pairs_task(run, det_name, flat_files, gains, mask_files=(),
                    flat2_finder=find_flat2_bot,
                    linearity_spec_range=(1e4, 9e4), use_exptime=False,
                    bias_frame=None, mondiode_func=None, dark_frame=None):
    """Single sensor execution of the flat pairs task."""
    file_prefix = make_file_prefix(run, det_name)

    task = sensorTest.FlatPairTask()
    task.run(file_prefix, flat_files, mask_files, gains,
             linearity_spec_range=linearity_spec_range,
             use_exptime=use_exptime, flat2_finder=flat2_finder,
             bias_frame=bias_frame, mondiode_func=mondiode_func,
             linearity_correction=get_nlc_func(det_name),
             dark_frame=dark_frame)

    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)

    detresp_file = '%s_det_response.fits' % file_prefix
    siteUtils.make_png_file(plots.linearity,
                            '%s_linearity.png' % file_prefix,
                            detresp_file=detresp_file, max_dev=0.03,
                            use_exptime=use_exptime,
                            Ne_bounds=linearity_spec_range)
    siteUtils.make_png_file(plots.linearity_resids,
                            '%s_linearity_resids.png' % file_prefix,
                            detresp_file=detresp_file, max_dev=0.03,
                            Ne_bounds=linearity_spec_range,
                            use_exptime=use_exptime)


def nonlinearity_task(run, det_name, detresp_file, outfile):
    """Single sensor execution of nonlinearity task."""
    file_prefix = make_file_prefix(run, det_name)
    task = sensorTest.NonlinearityTask()
    try:
        task.run(file_prefix, detresp_file, outputfile=outfile)
    except Exception:
        print(f'NonlinearityTask.run failed for {file_prefix}')


def ptc_task(run, det_name, flat_files, gains, mask_files=(),
             flat2_finder=find_flat2_bot, bias_frame=None):
    """Single sensor execution of the PTC task."""
    file_prefix = make_file_prefix(run, det_name)

    task = sensorTest.PtcTask()
    task.run(file_prefix, flat_files, mask_files, gains,
             flat2_finder=flat2_finder, bias_frame=bias_frame,
             linearity_correction=get_nlc_func(det_name))

    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)
    siteUtils.make_png_file(plots.ptcs,
                            '%s_ptcs.png' % file_prefix,
                            ptc_file='%s_ptc.fits' % file_prefix)


def bf_task(run, det_name, flat_files, gains, mask_files=(),
            flat2_finder=None, bias_frame=None):
    """Single sensor execution of the brighter-fatter task."""
    file_prefix = make_file_prefix(run, det_name)

    task = sensorTest.BFTask()
    task.run(file_prefix, flat_files, mask_files=mask_files,
             flat2_finder=flat2_finder, bias_frame=bias_frame,
             linearity_correction=get_nlc_func(det_name), gains=gains)

    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)
    siteUtils.make_png_file(plots.bf_curves,
                            '%s_brighter-fatter.png' % file_prefix,
                            bf_file='%s_bf.fits' % file_prefix)

def qe_jh_task(det_name):
    """JH version of single sensor execution of the QE task."""
    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    lambda_files = siteUtils.dependency_glob(glob_pattern('qe', det_name),
                                             acq_jobname=acq_jobname)
    if not lambda_files:
        print("qe_task: QE scan files not found for detector", det_name)
        return None

    pd_ratio_file = eotestUtils.getPhotodiodeRatioFile()
    if pd_ratio_file is None:
        message = ("The BOT photodiode ratio file is " +
                   "not given in config/%s/eotest_calibrations.cfg."
                   % siteUtils.getSiteName())
        raise RuntimeError(message)

#    correction_image = eotestUtils.getIlluminationNonUniformityImage()
#    if correction_image is None:
#        print()
#        print("WARNING: The correction image file is not given in")
#        print("config/%s/eotest_calibrations.cfg." % siteUtils.getSiteName())
#        print("No correction for non-uniform illumination will be applied.")
#        print()
#        sys.stdout.flush()
    mask_files = get_mask_files(det_name)
    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)
    bias_frame = bias_filename(run, det_name)

    return qe_task(run, det_name, lambda_files, pd_ratio_file, gains,
                   mask_files=mask_files, bias_frame=bias_frame,
                   mondiode_func=mondiode_value)

def qe_task(run, det_name, lambda_files, pd_ratio_file, gains,
            mask_files=(), correction_image=None, temp_set_point=-100,
            bias_frame=None, mondiode_func=None):
    """Single sensor execution of the QE task."""
    file_prefix = make_file_prefix(run, det_name)

    task = sensorTest.QeTask()
    task.config.temp_set_point = temp_set_point
    task.run(file_prefix, lambda_files, pd_ratio_file, mask_files, gains,
             correction_image=correction_image, bias_frame=bias_frame,
             mondiode_func=mondiode_func)

    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)

    siteUtils.make_png_file(plots.qe,
                            '%s_qe.png' % file_prefix,
                            qe_file='%s_QE.fits' % file_prefix)

    try:
        plots.flat_fields(os.path.dirname(lambda_files[0]),
                          annotation='e-/pixel, gain-corrected, bias-subtracted')
    except Exception as eobj:
        print("Exception raised while creating flat fields:")
        print(str(eobj))


def tearing_task(run, det_name, flat_files, bias_frame=None):
    """Single sensor execution of the tearing task."""
    file_prefix = make_file_prefix(run, det_name)

    tearing_found, _ = tearing_detection(flat_files, bias_frame=bias_frame)
    tearing_stats = [('BOT_EO_acq', 'N/A', det_name, len(tearing_found))]

    with open('%s_tearing_stats.pkl' % file_prefix, 'wb') as output:
        pickle.dump(tearing_stats, output)


def overscan_task(run, det_name, flat_files, gains, bias_frame=None):
    """Single sensor execution of the overscan task."""
    file_prefix = make_file_prefix(run, det_name)
    task = sensorTest.OverscanTask()
    task.run(file_prefix, flat_files, gains, bias_frame=bias_frame)

    overscan_file = f'{file_prefix}_overscan_results.fits'
    plots = OverscanTestPlots(file_prefix, overscan_file=overscan_file)

    plot_types = '''serial_eper_low serial_eper_high serial_cti
                    serial_overscan_signal serial_overscan_sum
                    parallel_eper_low parallel_eper_high parallel_cti
                    parallel_overscan_signal parallel_overscan_sum'''.split()

    plot_funcs = {_: getattr(plots, f'{_}_curves') for _ in plot_types}

    for plot_type, plot_func in plot_funcs.items():
        siteUtils.make_png_file(plot_func, f'{file_prefix}_{plot_type}.png')
        plt.close()


def persistence_task(run, det_name, bias_files, superbias_frame, mask_files):
    """Single sensor execution of the persistence analysis."""
    file_prefix = make_file_prefix(run, det_name)
    data = defaultdict(list)
    ccd = None
    for bias_file in bias_files:
        ccd = sensorTest.MaskedCCD(bias_file, mask_files=mask_files,
                                   bias_frame=superbias_frame)
        tseqnum = ccd.md.get('TSEQNUM')
        for amp in ccd:
            amp_image = ccd.unbiased_and_trimmed_image(amp)
            stats = afwMath.makeStatistics(amp_image,
                                           afwMath.MEAN | afwMath.STDEV,
                                           ccd.stat_ctrl)
            data['tseqnum'].append(tseqnum)
            data['amp'].append(amp)
            data['mean_signal'].append(stats.getValue(afwMath.MEAN))
            data['stdev'].append(stats.getValue(afwMath.STDEV))
    df = pd.DataFrame(data=data)
    outfile = f'{file_prefix}_persistence_data.pkl'
    df.to_pickle(outfile)
    fig = plt.figure()
    for amp in ccd:
        my_df = df.query(f'amp == {amp}')
        plt.scatter(my_df['tseqnum'], my_df['mean_signal'], s=2, label=f'{amp}')
    xmax = 1.2*(np.max(df['tseqnum']) - np.min(df['tseqnum'])) \
           + np.min(df['tseqnum'])
    axis = plt.axis()
    plt.xlim(axis[0], xmax)
    plt.legend(fontsize='x-small')
    plt.xlabel('test sequence number')
    plt.ylabel('mean residual signal (ADU)')
    plt.title(f'{file_prefix} persistence test')
    plt.savefig(f'{file_prefix}_persistence.png')
    plt.close()


def get_raft_files_by_slot(raft_name, file_suffix, jobname=None):
    """Return a dictionary of raft filenames, keyed by slot_name."""
    run = '*'
    template = '{}_{}_{}_' + file_suffix
    raft_files = dict()
    for slot_name in camera_info.get_slot_names():
        pattern = template.format(raft_name, slot_name, run)
        filenames = glob.glob(pattern)
        if filenames:
            raft_files[slot_name] = filenames[0]
        else:
            filenames = siteUtils.dependency_glob(pattern, jobname=jobname)
            if filenames:
                raft_files[slot_name] = filenames[0]
    if not raft_files:
        raise FileNotFoundError("no files found for raft %s with suffix %s"
                                % (raft_name, file_suffix))
    return raft_files


def repackage_summary_files():
    """
    Repackage summary.lims files from prior jobs as eotest results
    files.
    """
    run = siteUtils.getRunNumber()
    summary_files = siteUtils.dependency_glob('summary.lims')
    for det_name in camera_info.get_det_names():
        file_prefix = make_file_prefix(run, det_name)
        raft, slot = det_name.split('_')
        repackager = eotestUtils.JsonRepackager()
        repackager.process_files(summary_files, slot=slot, raft=raft)
        repackager.write('{}_eotest_results.fits'.format(file_prefix))


def get_bot_eo_config(bot_eo_config_file=None):
    """Get the ConfigParser object from the BOT EO configuration file."""
    bot_eo_config_file = siteUtils.get_bot_eo_config_file(bot_eo_config_file)

    # Read in the analyses to be performed from the config file.
    cp = configparser.ConfigParser(allow_no_value=True,
                                   inline_comment_prefixes=("#", ))
    cp.optionxform = str   # allow for case-sensitive keys
    cp.read(bot_eo_config_file)
    return cp


def get_analysis_types(bot_eo_config_file=None):
    """"Get the analysis types to be performed from the BOT-level EO config."""
    cp = get_bot_eo_config(bot_eo_config_file)

    analysis_types = []
    if 'ANALYZE' not in cp:
        return analysis_types
    for analysis_type, _ in cp.items("ANALYZE"):
        analysis_types.append(analysis_type)

    return analysis_types


def get_nlc_func(det_name, bot_eo_config_file=None):
    """
    Return the nonlinearity correction function for the specified
    detector.
    """
    file_prefix = make_file_prefix('*', det_name)
    try:
        nlc_file = siteUtils.dependency_glob(f'{file_prefix}_nlc.fits',
                                             jobname='nonlinearity_BOT')[0]
    except IndexError as eobj:
        print(f'{file_prefix}_nlc.fits not found:', eobj)
        return None

    cp = get_bot_eo_config(bot_eo_config_file)
    if 'NLC_PARAMS' in cp:
        nlc_params = {k: siteUtils.cast(v) for k, v in cp['NLC_PARAMS'].items()}
    else:
        nlc_params = dict()

    try:
        nlc = sensorTest.NonlinearityCorrection.create_from_fits_file(
            nlc_file, **nlc_params)

    except ValueError as eobj:
        print(f'Error creating nlc from {nlc_file}:', eobj)
        nlc = None

    return nlc


def run_jh_tasks(*jh_tasks, device_names=None, processes=None, walltime=3600):
    """
    Run functions to execute tasks under the job harness in parallel.
    These functions should take a device name as its only argument, and
    the parallelization will take place over device_names.


    Parameters
    ----------
    jh_tasks: list-like container of functions
        These functions are serialized and dispatched to workers on
        remote nodes, so all dependencies should be imported in the
        bodies of the functions.
    device_names: list-like container of device names [None]
        List of sensors or rafts on which to operate.  If None, then
        the installed sensors in the focal plane is used.
    processes: int [None]
        Number of processes to run in parallel. If None, then all
        available processes can be potentially used.
    walltime: float [3600]
        Walltime in seconds for python app execution.  If the python app
        does not return within walltime, a parsl.app.errors.AppTimeout
        exception will be thrown.

    Raises
    ------
    parsl.app.errors.AppTimeout


    Notes
    -----
    Because the number of jh_task functions can vary, the keyword arguments
    should reference the keywords explicitly, i.e., one cannot rely on
    keyword position to pass those values.
    """
    if device_names is None:
        device_names = camera_info.get_det_names()

    # Restrict to installed rafts or sensors.  This function call
    # also writes the camera_info cache file for the eT db query.
    installed_rafts = camera_info.get_installed_raft_names()

    # Check if rafts are over-ridden in the lcatr.cfg file.
    override_rafts = os.environ.get('LCATR_RAFTS', None)
    if override_rafts is not None:
        installed_rafts = override_rafts.split('_')

    device_names = [_ for _ in device_names if _[:3] in installed_rafts]

    cwd = os.path.abspath('.')

    # Query eT database for file paths from a previous run, if
    # specified, and store in a pickle file.
    hj_fp_server = siteUtils.HarnessedJobFilePaths()

    # Query for file paths for other analysis runs, if specified in
    # the bot_eo_config_file.
    for analysis_type in ('badpixel', 'bias', 'dark', 'linearity',
                          'nonlinearity'):
        hj_fp_server.query_file_paths(
            siteUtils.get_analysis_run(analysis_type))

    hj_fp_server_file = 'hj_fp_server.pkl'
    with open(hj_fp_server_file, 'wb') as output:
        pickle.dump(hj_fp_server, output)

    # Create a GetAmplifierGains object in order to query the eT
    # database for gain results from previous runs and write a pickle
    # file that can be loaded locally from disk by the various jh
    # tasks being run in parallel to avoid eT db access contention.
    GetAmplifierGains()

    for jh_task in jh_tasks:
        # Add 30 second sleep before launching jh_task processes in
        # parallel to allow for parsl process_pool_workers from the
        # previous set of jh_task processes to finish.
        time.sleep(30)
        run_device_analysis_pool(jh_task, device_names,
                                 processes=processes, cwd=cwd,
                                 walltime=walltime)


def run_python_task_or_cl_script(python_task, cl_script, device_names=None,
                                 processes=None, walltime=3600):
    """
    If we are running a traveler for the Cryostat, use the
    ssh_dispatcher to run jobs in parallel on the diagnostic cluster,
    otherwsie we are running at TS8, so use multiprocessing on the
    current node.
    """
    # Copy the bot_eo_acq_cfg file to the current working directory
    # so it can be persisted in the eT and DataCatalog.
    acq_config = siteUtils.get_job_acq_configs()
    bot_eo_acq_cfg = acq_config['bot_eo_acq_cfg']
    shutil.copy(bot_eo_acq_cfg, '.')

    if (os.environ.get('LCATR_USE_PARSL', False) == 'True'
        or siteUtils.getUnitType() == 'LCA-10134_Cryostat'):
        # Run command-line verions using parsl or ssh_dispatcher.
        run_jh_tasks(cl_script, device_names=device_names,
                     processes=processes, walltime=walltime)
    else:
        # Run python version using multiprocessing directly.
        run_jh_tasks(python_task, device_names=device_names,
                     processes=processes, walltime=walltime)
