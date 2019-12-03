#!/usr/bin/env ipython
"""
Script for BOT flat gain stability analysis.
"""
def flat_gain_stability_jh_task(det_name):
    """JH version of single sensor execution of the flat pairs task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        bias_filename, flat_gain_stability_task,\
        get_mask_files, medianed_dark_frame

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    flat_files = siteUtils.dependency_glob(glob_pattern('tearing', det_name),
                                           acq_jobname=acq_jobname)
    if not flat_files:
        print("flat_gain_stability_task: Flat pairs files not found for",
              det_name)
        return None

    mask_files = get_mask_files(det_name)
    bias_frame = bias_filename(run, det_name)
    dark_frame = medianed_dark_frame(det_name)

    return flat_gain_stability_task(run, det_name, flat_files,
                                    mask_files=mask_files,
                                    bias_frame=bias_frame,
                                    dark_frame=dark_frame)


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    flat_gain_stability_jh_task(det_name)
