#!/usr/bin/env python
"""
Producer script for BOT tearing analysis.
"""
def tearing_jh_task(det_name):
    """JH version of single sensor execution of the tearing task."""
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        bias_filename, tearing_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    flat_files = siteUtils.dependency_glob(glob_pattern('tearing', det_name),
                                           acq_jobname=acq_jobname)
    if not flat_files:
        print("tearing_task: Flat files not found for detector", det_name)
        return None
    bias_frame = bias_filename(file_prefix)
    return tearing_task(run, det_name, flat_files, bias_frame=bias_frame)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'tearing' in get_analysis_types():
        run_jh_tasks(tearing_jh_task)
