#!/usr/bin/env python
"""
Producer script for BOT flat pairs (linearity and full-well) analysis.
"""
def flat_pairs_jh_task(det_name):
    """JH version of single sensor execution of the flat pairs task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, flat_pairs_task, mondiode_value

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    flat_files \
        = siteUtils.dependency_glob(glob_pattern('flat_pairs', det_name),
                                    acq_jobname=acq_jobname)
    if not flat_files:
        print("flat_pairs_task: Flat pairs files not found for detector",
              det_name)
        return None

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))
    try:
        eotest_results_file \
            = siteUtils.dependency_glob('{}_eotest_results.fits'.format(file_prefix),
                                        jobname='fe55_analysis_BOT')[0]
    except IndexError:
        print("flat_pairs_jh_task: Fe55 eotest results file not found for ",
              file_prefix, ".  Using unit gains.")
        eotest_results_file = None
    gains = get_amplifier_gains(eotest_results_file)
    bias_frame = bias_filename(file_prefix)

    return flat_pairs_task(run, det_name, flat_files, gains,
                           mask_files=mask_files, bias_frame=bias_frame,
                           mondiode_func=mondiode_value)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'linearity' in get_analysis_types():
        run_jh_tasks(flat_pairs_jh_task)
