#!/usr/bin/env python
"""
Producer script for BOT brighter-fatter analysis.
"""
def bf_jh_task(det_name):
    """JH version of single sensor execution of the brighter-fatter task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        bias_filename, bf_task, find_flat2_bot

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    flat_files \
        = siteUtils.dependency_glob(glob_pattern('brighter_fatter', det_name),
                                    acq_jobname=acq_jobname)

    if not flat_files:
        print("bf_jh_task: Flat pairs files not found for detector", det_name)
        return None

    flat_files = [_ for _ in flat_files if 'flat1' in _]

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))
    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    bias_frame = bias_filename(file_prefix)

    return bf_task(run, det_name, flat_files, mask_files=mask_files,
                   flat2_finder=find_flat2_bot, bias_frame=bias_frame)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'brighterfatter' in get_analysis_types():
        run_jh_tasks(bf_jh_task)
