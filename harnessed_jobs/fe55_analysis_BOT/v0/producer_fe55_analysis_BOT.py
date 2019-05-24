#!/usr/bin/env python
"""
Producer script for BOT Fe55 analysis.
"""
def fe55_jh_task(det_name):
    "JH version of single sensor execution of the Fe55 analysis task."
    import os
    import siteUtils
    from bot_eo_analyses import glob_pattern, fe55_task

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    nbias = os.environ.get('LCATR_NUM_BIAS_FRAMES', 25)

    fe55_files = siteUtils.dependency_glob(glob_pattern('fe55', det_name),
                                           acq_jobname=acq_jobname)
    bias_files = siteUtils.dependency_glob(glob_pattern('fe55_bias', det_name),
                                           acq_jobname=acq_jobname)[:nbias]
    if not fe55_files or not bias_files:
        print("fe55_task: Needed data files missing for detector", det_name)
        return None
    return fe55_task(run, det_name, fe55_files, bias_files)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'gain' in get_analysis_types():
        run_jh_tasks(fe55_jh_task)
