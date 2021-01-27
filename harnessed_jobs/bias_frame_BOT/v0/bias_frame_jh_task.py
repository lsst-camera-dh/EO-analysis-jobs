#!/usr/bin/env ipython
"""
Command-line script for detector-level bias frame job harness task.
"""
def bias_frame_jh_task(det_name):
    """JH version of the bias_frame_task."""
    import os
    import siteUtils
    import json
    from bot_eo_analyses import glob_pattern, bias_frame_task, \
        bias_stability_task

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    bias_files \
        = siteUtils.dependency_glob(glob_pattern('bias_frame', det_name),
                                    acq_jobname=acq_jobname,
                                    description='Bias frames:')
    if not bias_files:
        print("bias_frame_task: Needed data files are missing for detector",
              det_name)
        return None

    bias_stability_files \
        = siteUtils.dependency_glob(glob_pattern('bias_stability', det_name),
                                    acq_jobname=acq_jobname,
                                    description='Bias stability frames:')
    if not bias_stability_files:
        print("bias_stability_task: Needed data files are missing for detector",
              det_name)
        return None

    bias_frame, pca_files = bias_frame_task(run, det_name, bias_files[1:])
    bias_stability_files = sorted(bias_stability_files)

    if not os.environ.get('LCATR_USE_PCA_BIAS_FIT', False) == 'True':
        pca_files = None
    print("pca_files:", pca_files)
    bias_stability_task(run, det_name, bias_stability_files,
                        pca_files=pca_files)

    return bias_frame


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    bias_frame_jh_task(det_name)
