#!/usr/bin/env ipython
"""
Script for BOT dark current analysis.
"""
def dark_current_jh_task(det_name):
    """JH version of single sensor execution of the dark current task."""
    from collections import defaultdict
    from astropy.io import fits
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, dark_current_task,\
        plot_ccd_total_noise, get_mask_files
    from bot_data_handling import most_common_dark_files

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    dark_files \
        = siteUtils.dependency_glob(glob_pattern('dark_current', det_name),
                                    acq_jobname=acq_jobname,
                                    description="Dark current frames:")
    if not dark_files:
        print("dark_current_task: No dark files found for detector", det_name)
        return None

    dark_files_linear_fit = list(dark_files)
    dark_files = most_common_dark_files(dark_files)
    if len(dark_files_linear_fit) == len(dark_files):
        # These data only have one integration time, so skip linear
        # fit of dark current signal vs integration time.
        dark_files_linear_fit = None

    mask_files = get_mask_files(det_name)
    eotest_results_file \
        = siteUtils.dependency_glob('{}_eotest_results.fits'.format(file_prefix),
                                    jobname='read_noise_BOT')[0]
    gains = get_amplifier_gains('{}_eotest_results.fits'.format(file_prefix))
    bias_frame = bias_filename(run, det_name)

    dark_curr_pixels, dark95s \
        = dark_current_task(run, det_name, dark_files, gains,
                            mask_files=mask_files, bias_frame=bias_frame,
                            dark_files_linear_fit=dark_files_linear_fit)
    plot_ccd_total_noise(run, det_name, dark_curr_pixels, dark95s,
                         eotest_results_file)
    return dark_curr_pixels, dark95s


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    dark_current_jh_task(det_name)
