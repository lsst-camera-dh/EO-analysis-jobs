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

    png_files = ['%s_gains.png' % file_prefix]
    siteUtils.make_png_file(plots.gains, png_files[-1])

    png_files.append('%s_mean_bias.png' % file_prefix)
    siteUtils.make_png_file(sensorTest.plot_flat, png_files[-1], mean_bias_file,
                            title='%s, mean bias frame' % det_name,
                            annotation='ADU/pixel, overscan-subtracted')

    fe55_file = glob.glob('%s_psf_results*.fits' % file_prefix)[0]
    png_files.append('%s_fe55_dists.png' % file_prefix)
    siteUtils.make_png_file(plots.fe55_dists, png_files[-1],
                            fe55_file=fe55_file)

    png_files.append('%s_psf_dists.png' % file_prefix)
    siteUtils.make_png_file(plots.psf_dists, png_files[-1],
                            fe55_file=fe55_file)

    with open('fe55_task_png_files.txt', 'w') as output:
        for item in png_files:
            output.write('{}\n'.format(item))


def get_amplifier_gains(eotest_results_file):
    """Extract the gains for each amp in an eotest_results file."""
    data = sensorTest.EOTestResults(eotest_results_file)
    amps = data['AMP']
    gains = data['GAIN']
    return {amp: gain for amp, gain in zip(amps, gains)}


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

    title = '%s, medianed dark for bright defects analysis' % file_prefix
    annotation='e-/pixel, gain-corrected, bias-subtracted'
    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_medianed_dark.png' % file_prefix,
                            '%s_median_dark_bp.fits' % file_prefix,
                            title=title, annotation=annotation)


def dark_defects_task(det_name):
    """Single sensor execution of the dark defects task."""
    run = siteUtils.getRunNumber()
    file_prefix = '%s_%s' % (run, det_name)
    pattern = 'sflat_500_flat_H/*_{}.fits'.format(det_name)
    sflat_files = siteUtils.dependency_glob(pattern)
    if not sflat_files:
        print("No high flux superflat files found for run", run)
        return
    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    task = sensorTest.DarkPixelsTask()
    task.run(det_name, sflat_files, mask_files)

    title = '%s, superflat for dark defects analysis' % file_prefix
    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_superflat_dark_defects.png' % file_prefix,
                            '%s_median_sflat.fits' % file_prefix,
                            title=title, annotation='ADU/pixel')


def traps_task(det_name):
    """Single sensor execution of the traps analysis task."""
    run = siteUtils.getRunNumber()
    file_prefix = '%s_%s' % (run, det_name)
    pattern = 'trap_ppump_000/*_{}.fits'.format(det_name)
    trap_files = siteUtils.dependency_glob(pattern))
    if not trap_files:
        print("No pocket pumping file found for run", run)
        return
    trap_file = trap_files[0]
    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    # Omit rolloff defects mask since a trap in the rolloff edge region can
    # affect the entire column.
    mask_files \
        = [item for item in mask_files if item.find('edge_rolloff') == -1]

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    task = sensorTest.TrapTask()
    task.run(det_name, trap_file, mask_files, gains)


def dark_current_task(det_name):
    """Single sensor execution of the dark current task."""
    run = siteUtils.getRunNumber()
    file_prefix = '%s_%s' % (run, det_name)

    pattern = 'dark_dark_*/*_{}.fits'.format(det_name)
    dark_files = siteUtils.dependency_glob(pattern))
    if not dark_files:
        print("No dark files found for run:", run)
        return

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(det_name)))

    eotest_results_file = '{}_eotest_results.fits'.format(det_name)
    gains = get_amplifier_gains(eotest_results_file)

    task = sensorTest.DarkCurrentTask()
    task.config.temp_set_point = -100.
    dark_curr_pixels, dark95s \
        = task.run(det_name, dark_files, mask_files, gains)

    eo_results = sensorTest.EOTestResults(eotest_results_file)
    read_noise = dict(pair for pair in zip(eo_results['AMP'],
                                           eo_results['TOTAL_NOISE']))

    siteUtils.make_png_file(sensorTest.total_noise_histograms,
                            '%s_total_noise_hists.png' % file_prefix,
                            dark_curr_pixels, read_noise, dark95s,
                            exptime=16, title=det_name)

    plots = sensorTest.EOTestPlots(det_name, results_file=eotest_results_file)
    siteUtils.make_png_file(plots.total_noise, '%s_noise.png' % file_prefix,
                            dark95s=dark95s)


if __name__ == '__main__':
    det_names = camera_info.get_det_names()
    raft_names = camera_info.get_raft_names()

    processes = None
    run_sensor_analysis_pool(fe55_task, det_names, processes=processes)
    run_sensor_analysis_pool(read_noise_task, det_names, processes=processes)
    for raft_name in raft_names:
        raft_noise_correlations(raft_name)
    run_sensor_analysis_pool(bright_defects_task, det_names,
                             processes=processes)
    run_sensor_analysis_pool(dark_defects_task, det_names, processes=processes)
    run_sensor_analysis_pool(traps_task, det_names, processes=processes)
    run_sensor_analysis_pool(dark_current_task, det_names, processes=processes)
