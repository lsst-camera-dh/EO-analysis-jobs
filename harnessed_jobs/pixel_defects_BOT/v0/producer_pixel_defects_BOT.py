#!/usr/bin/env python
"""
Producer script for BOT bright and dark pixel analyses.
"""
def bright_defects_jh_task(det_name):
    """JH version of single sensor bright pixels task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, bright_defects_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    dark_files \
        = siteUtils.dependency_glob(glob_pattern('bright_defects', det_name),
                                    acq_jobname=acq_jobname)
    if not dark_files:
        print("bright_defects_task: Needed data files missing for detector",
              det_name)
        return None

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)
    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    bias_frame = bias_filename(file_prefix)

    return bright_defects_task(run, det_name, dark_files, gains,
                               mask_files=mask_files, bias_frame=bias_frame)


def dark_defects_jh_task(det_name):
    """JH version of single sensor execution of the dark defects task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, dark_defects_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    sflat_files \
        = siteUtils.dependency_glob(glob_pattern('dark_defects', det_name),
                                    acq_jobname=acq_jobname)
    if not sflat_files:
        print("dark_defects_task: No high flux superflat files found for",
              det_name)
        return None

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))
    bias_frame = bias_filename(file_prefix)

    return dark_defects_task(run, det_name, sflat_files, mask_files=mask_files,
                             bias_frame=bias_frame)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'badpixel' in get_analysis_types():
        run_jh_tasks(bright_defects_jh_task, dark_defects_jh_task)
