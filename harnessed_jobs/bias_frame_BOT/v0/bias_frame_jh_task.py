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
        bias_stability_task, pca_corrected_superbias

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    bias_files \
        = siteUtils.dependency_glob(glob_pattern('bias_frame', det_name),
                                    acq_jobname=acq_jobname,
                                    description='Bias frames:')
    # Skip the first 5 bias files to avoid transient features in the
    # imaging array from after just starting the run.
    bias_files = sorted(bias_files, key=os.path.basename)[5:]
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

    superbias_file, bias_model_components \
        = bias_frame_task(run, det_name, bias_files)
    bias_stability_files = sorted(bias_stability_files)

    if not os.environ.get('LCATR_USE_PCA_BIAS_FIT', "True") == 'True':
        bias_model_components = None
    if siteUtils.get_analysis_run('bias') == 'rowcol':
        bias_model_components = 'rowcol', superbias_file
    print("bias_model_components:", bias_model_components)
    bias_stability_task(run, det_name, bias_stability_files,
                        bias_model_components=bias_model_components)

    if (bias_model_components is not None and
        bias_model_components[0] != 'rowcol'):
        pca_corrected_superbias(run, det_name, bias_files,
                                bias_model_components)

    return superbias_file


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    bias_frame_jh_task(det_name)
