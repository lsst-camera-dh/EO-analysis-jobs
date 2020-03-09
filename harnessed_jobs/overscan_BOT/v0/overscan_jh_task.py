#!/usr/bin/env ipython
"""
Harnessed job script for BOT overscan analysis.
"""
def overscan_jh_task(det_name):
    """JH version of single sensor execution of the overscan task."""
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, overscan_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    flat_files = siteUtils.dependency_glob(glob_pattern('overscan', det_name),
                                           acq_jobname=acq_jobname)
    if not flat_files:
        print("overscan_task: Flat pairs files not found for detector",
              det_name)
        return None

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)
    bias_frame = bias_filename(run, det_name)

    return overscan_task(run, det_name, flat_files, gains,
                         bias_frame=bias_frame)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    overscan_jh_task(det_name)
