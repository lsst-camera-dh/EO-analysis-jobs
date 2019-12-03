#!/usr/bin/env ipython
"""
Producer script for BOT gain stability analysis.  The fe55_jh_task.py
script needs to be executed before this script.
"""
def gain_stability_jh_task(det_name):
    "JH version of the gain stability analysis task."
    import os
    import siteUtils
    from bot_eo_analyses import glob_pattern, gain_stability_task

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    fe55_files = siteUtils.dependency_glob(glob_pattern('fe55', det_name),
                                           acq_jobname=acq_jobname)
    if not fe55_files:
        print("fe55_task: Needed data files missing for detector", det_name)
        return None
    return gain_stability_task(run, det_name, fe55_files)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    gain_stability_jh_task(det_name)
