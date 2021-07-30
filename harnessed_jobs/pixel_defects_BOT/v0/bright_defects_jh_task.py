#!/usr/bin/env ipython
"""
Script for BOT bright pixel analysis.
"""
def bright_defects_jh_task(det_name):
    """JH version of single sensor bright pixels task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, bright_defects_task
    from bot_data_handling import most_common_dark_files

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    dark_files \
        = siteUtils.dependency_glob(glob_pattern('bright_defects', det_name),
                                    acq_jobname=acq_jobname)
    dark_files = most_common_dark_files(dark_files)
    if not dark_files:
        print("bright_defects_task: Needed data files missing for detector",
              det_name)
        return None

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)
    mask_files = sorted(glob.glob(f'{file_prefix}*mask*.fits'))
    bias_frame = bias_filename(run, det_name)

    return bright_defects_task(run, det_name, dark_files, gains,
                               mask_files=mask_files, bias_frame=bias_frame)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    bright_defects_jh_task(det_name)
