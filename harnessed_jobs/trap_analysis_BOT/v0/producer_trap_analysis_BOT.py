#!/usr/bin/env python
"""
Producer script for BOT trap analysis.
"""
def traps_jh_task(det_name):
    """JH version of single sensor execution of the traps analysis task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, traps_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    trap_files = siteUtils.dependency_glob(glob_pattern('traps', det_name),
                                           acq_jobname=acq_jobname)
    if not trap_files:
        print("traps_task: No pocket pumping file found for detector",
              det_name)
        return None
    trap_file = trap_files[0]

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    # Omit rolloff defects mask since a trap in the rolloff edge region can
    # affect the entire column.
    mask_files \
        = [item for item in mask_files if item.find('edge_rolloff') == -1]

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    bias_frame = bias_filename(file_prefix)

    return traps_task(run, det_name, trap_file, gains, mask_files=mask_files,
                      bias_frame=bias_frame)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'traps' in get_analysis_types():
        run_jh_tasks(traps_jh_task)
