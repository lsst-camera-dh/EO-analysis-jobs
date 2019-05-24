#!/usr/bin/env python
"""
Producer script for BOT bias frame generation.
"""
def bias_frame_jh_task(det_name):
    """JH version of the bias_frame_task."""
    import os
    import siteUtils
    from bot_eo_analyses import glob_pattern, bias_frame_task

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    nbias = os.environ.get('LCATR_NUM_BIAS_FRAMES', 10)
    bias_files \
        = siteUtils.dependency_glob(glob_pattern('bias_frame', det_name),
                                    acq_jobname=acq_jobname)[:nbias]
    if not bias_files:
        print("bias_frame_task: Needed data files are missing for detector",
              det_name)
        return None
    return bias_frame_task(run, det_name, bias_files)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'bias' in get_analysis_types():
        run_jh_tasks(bias_frame_jh_task)
