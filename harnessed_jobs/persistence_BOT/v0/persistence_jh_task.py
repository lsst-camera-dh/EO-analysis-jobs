#!/usr/bin/env ipython
"""
Command-line script for detector-level bias frame job harness task.
"""
def persistence_jh_task(det_name):
    """JH version of the persistence_task."""
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern, \
        bias_frame_task, get_mask_files, persistence_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)

    acq_jobname = siteUtils.getProcessName('BOT_acq')
    bias_files \
        = siteUtils.dependency_glob(glob_pattern('persistence_bias', det_name),
                                    acq_jobname=acq_jobname,
                                    description='Persistence bias frames:')
    if not bias_files:
        print("persistence_task: Needed data files are missing for detector",
              det_name)
        return None

    # We should parse this number from the bot_eo_acq_config file:
    num_post_exp_bias_frames = 5

    # Make a superbias frame using the pre-exposure persistence bias
    # files, skipping the first exposure.
    superbias_frame = f'{file_prefix}_persistence_superbias.fits'
    bias_frame_task(run, det_name, bias_files[1:-num_post_exp_bias_frames],
                    bias_frame=superbias_frame)

    return persistence_task(run, det_name, bias_files[-num_post_bias_frames:],
                            superbias_frame, get_mask_files(det_name))

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    persistence_jh_task(det_name)
