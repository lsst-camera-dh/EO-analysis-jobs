#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
from __future__ import print_function
import glob
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import siteUtils
from correlated_noise import correlated_noise, raft_level_oscan_correlations
from camera_components import camera_info
from multiprocessor_execution import run_sensor_analysis_pool


def get_amplifier_gains(eotest_results_file):
    """Extarct the gains for each amp in an eotest_results file."""
    data = sensorTest.EOTestResults(eotest_results_file)
    amps = data['AMP']
    gains = data['GAIN']
    return {amp: gain for amp, gain in zip(amps, gains)}


def bright_defects_task(det_name):
    "Single sensor execution of the bright pixels task."
    run = siteUtils.getRunNumber()
    file_prefix = '%s_%s' % (run, det_name)
    dark_files \
        = siteUtils.dependency_glob('dark_dark_*/*_{}.fits'.format(det_name))
    if not dark_files:
        print("Needed data files missing for run {}, detector {}."
              .format(run, det_name))
        return

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)
    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    task = sensorTest.BrightPixelsTask()
    task.config.temp_set_point = -100.
    task.run(file_prefix, dark_files, mask_files, gains)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_medianed_dark.png' % file_prefix,
                            '%s_median_dark_bp.fits' % file_prefix,
                            title='%s, medianed dark for bright defects analysis' % file_prefix,
                            annotation='e-/pixel, gain-corrected, bias-subtracted')


def read_noise_task(det_name):
    """Run the single sesnor read noise task."""
    run = siteUtils.getRunNumber()
    file_prefix = '%s_%s' % (run, det_name)
    bias_files \
        = siteUtils.dependency_glob('fe55_fe55_*/*_{}.fits'.format(det_name))
    if not bias_files:
        print("Needed data files are missing for run {}, detector {}."
              .format(run, det_name))
        return
    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    system_noise = None

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    task = sensorTest.ReadNoiseTask()
    task.config.temp_set_point = -100.
    task.run(file_prefix, bias_files, gains, system_noise=system_noise,
             mask_files=mask_files, use_overscan=True)

    # Compute amp-amp correlated noise.
    _, corr_fig, _ = correlated_noise(bias_files, target=0,
                                      make_plots=True, title=file_prefix)
    plt.figure(corr_fig.number)
    plt.savefig('%s_correlated_noise.png' % file_prefix)


def raft_noise_correlations(raft_name):
    """This is the raft-level noise-correlation analysis."""
    run = siteUtils.getRunNumber()
    file_prefix = '{}_{}'.format(run, raft_name)
    bias_files = dict()
    slot_names = camera_info.get_slot_names()
    for slot_name in slot_names:
        det_name = '{}_{}'.format(raft_name, slot_name)
        my_files \
            = siteUtils.dependency_glob('fe55_bias_*/*_{}.fits'.format(det_name))
        if not my_files:
            print("Missing bias files for raft {}. Skipping.".format(raft_name))
            return
        bias_files[slot_name] = my_files[0]
    title = 'Overscan correlations, Run {}, {}'.format(raft_name, run)
    raft_level_oscan_correlations(bias_files, title=title)
    plt.savefig('{}_overscan_correlations.png'.format(file_prefix))


def fe55_task(det_name):
    "Single sensor execution of the Fe55 analysis task."
    fe55_files \
        = siteUtils.dependency_glob('fe55_fe55_*/*_{}.fits'.format(det_name))
    bias_files \
        = siteUtils.dependency_glob('fe55_bias_*/*_{}.fits'.format(det_name))
    file_prefix = '{}_{}'.format(siteUtils.getRunNumber(), det_name)

    if not fe55_files or not bias_files:
        print("Needed data files missing for run {}, detector {}."
              .format(siteUtils.getRunNumber(), det_name))
        return

    mean_bias_file = '{}_mean_bias.fits'.format(file_prefix)
    imutils.fits_mean_file(bias_files, mean_bias_file)

    pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=det_name)

    siteUtils.make_png_file(pixel_stats.pixel_hists,
                            '%s_fe55_p3_p5_hists.png' % file_prefix,
                            pix0='p3', pix1='p5')

    siteUtils.make_png_file(pixel_stats.pixel_diff_profile,
                            '%s_fe55_p3_p5_profiles.png' % file_prefix,
                            pixel_coord='x', pix0='p3', pix1='p5')

    siteUtils.make_png_file(pixel_stats.apflux_profile,
                            '%s_fe55_apflux_serial.png' % file_prefix)

    siteUtils.make_png_file(pixel_stats.apflux_profile,
                            '%s_fe55_apflux_parallel.png' % file_prefix,
                            pixel_coord='y')

    rolloff_mask_file = '%s_edge_rolloff_mask.fits' % file_prefix
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    task = sensorTest.Fe55Task()
    task.config.temp_set_point = -100.
    task.run(file_prefix, fe55_files, (rolloff_mask_file,), accuracy_req=0.01)

    # Fe55 gain and psf analysis results plots for the test report.
    results_file = '%s_eotest_results.fits' % file_prefix
    plots = sensorTest.EOTestPlots(file_prefix, results_file=results_file)

    siteUtils.make_png_file(plots.gains, '%s_gains.png' % file_prefix)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_mean_bias.png' % file_prefix,
                            mean_bias_file,
                            title='%s, mean bias frame' % det_name,
                            annotation='ADU/pixel, overscan-subtracted')

    fe55_file = glob.glob('%s_psf_results*.fits' % file_prefix)[0]
    siteUtils.make_png_file(plots.fe55_dists, '%s_fe55_dists.png' % file_prefix,
                            fe55_file=fe55_file)

    siteUtils.make_png_file(plots.psf_dists, '%s_psf_dists.png' % file_prefix,
                            fe55_file=fe55_file)

if __name__ == '__main__':
    det_names = camera_info.get_det_names()
    raft_names = camera_info.get_raft_names()

    run_sensor_analysis_pool(fe55_task, det_names, processes=None)
    run_sensor_analysis_pool(read_noise_task, det_names, processes=None)
    for raft_name in raft_names:
        raft_noise_correlations(raft_name)
    run_sensor_analysis_pool(bright_defects_task, det_names, processes=None)
