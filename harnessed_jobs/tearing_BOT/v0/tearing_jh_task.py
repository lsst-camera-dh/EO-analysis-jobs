#!/usr/bin/env python
"""
Script for BOT tearing analysis.
"""
def tearing_jh_task(det_name):
    """JH version of single sensor execution of the tearing task."""
    import siteUtils
    from bot_eo_analyses import glob_pattern, bias_filename, tearing_task

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    flat_files = siteUtils.dependency_glob(glob_pattern('tearing', det_name),
                                           acq_jobname=acq_jobname)
    if not flat_files:
        print("tearing_task: Flat files not found for detector", det_name)
        return None
    bias_frame = bias_filename(run, det_name)
    return tearing_task(run, det_name, flat_files, bias_frame=bias_frame)


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    tearing_jh_task(det_name)
