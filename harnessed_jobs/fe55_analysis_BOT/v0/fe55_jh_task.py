#!/usr/bin/env ipython
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

    fe55_files = siteUtils.dependency_glob(glob_pattern('fe55', det_name),
                                           acq_jobname=acq_jobname)
    if not fe55_files:
        print("fe55_task: Needed data files missing for detector", det_name)
        return None
    return fe55_task(run, det_name, fe55_files)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    fe55_jh_task(det_name)
