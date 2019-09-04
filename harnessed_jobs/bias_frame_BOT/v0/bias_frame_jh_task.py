#!/usr/bin/env python
"""
Command-line script for detector-level bias frame job harness task.
"""
def bias_frame_jh_task(det_name):
    """JH version of the bias_frame_task."""
    import os
    import siteUtils
    import json
    from bot_eo_analyses import glob_pattern, bias_frame_task
    if os.environ.get('LCATR_USE_PARSL', 'False') == 'True':
        suffix = 'with_parsl'
        import parsl
        print("parsl:", parsl)
    else:
        suffix = 'without_parsl'

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
    import sys
    det_name = sys.argv[1]
    bias_frame_jh_task(det_name)
