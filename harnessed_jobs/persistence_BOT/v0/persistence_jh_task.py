#!/usr/bin/env ipython
"""
Command-line script for detector-level bias frame job harness task.
"""
def persistence_jh_task(det_name):
    """JH version of the persistence_task."""
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern, \
        bias_frame_task, get_mask_files, get_bot_eo_config, persistence_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)

    acq_jobname = siteUtils.getProcessName('BOT_acq')
    bias_files \
        = siteUtils.dependency_glob(glob_pattern('persistence_bias', det_name),
                                    acq_jobname=acq_jobname,
                                    description='Persistence bias frames:')
    dark_files \
        = siteUtils.dependency_glob(glob_pattern('persistence_dark', det_name),
                                    acq_jobname=acq_jobname,
                                    description='Persistence dark frames:')
    if not bias_files or not dark_files:
        print("persistence_task: Needed data files are missing for detector",
              det_name)
        return None

    # Sort by test sequence number, i.e., by filenames.
    bias_files = sorted(bias_files)
    dark_files = sorted(dark_files)

    # Make a superbias frame using the pre-exposure persistence bias
    # files, skipping the first exposure.
    superbias_frame = f'{file_prefix}_persistence_superbias.fits'
    bias_frame_task(run, det_name, bias_files, bias_frame=superbias_frame)

    return persistence_task(run, det_name, dark_files, superbias_frame,
                            get_mask_files(det_name))


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    persistence_jh_task(det_name)
