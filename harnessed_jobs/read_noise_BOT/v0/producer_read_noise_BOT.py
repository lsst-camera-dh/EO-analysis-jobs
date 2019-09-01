#!/usr/bin/env python
"""
Producer script for BOT read noise analysis.
"""
def read_noise_jh_task(det_name):
    """JH version of the single sensor read noise task."""
    import os
    import glob
    import logging
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, read_noise_task

    logger = logging.getLogger('read_noise_jh_task')
    logger.setLevel(logging.INFO)

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    nbias = os.environ.get('LCATR_NUM_BIAS_FRAMES', 10)

    bias_files \
        = siteUtils.dependency_glob(glob_pattern('read_noise', det_name),
                                    acq_jobname=acq_jobname)[:nbias]
    if not bias_files:
        logger.info("read_noise_task: Needed data files are missing "
                    "for detector %s", det_name)
        return None
    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    return read_noise_task(run, det_name, bias_files, gains,
                           mask_files=mask_files)


def raft_jh_noise_correlations(raft_name):
    """JH version of raft-level noise-correlation analysis."""
    import os
    from collections import defaultdict
    import logging
    import siteUtils
    from camera_components import camera_info
    from bot_eo_analyses import raft_noise_correlations, glob_pattern

    logger = logging.getLogger('raft_jh_noise_correlations')
    logger.setLevel(logging.INFO)

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    bias_files \
        = siteUtils.dependency_glob(glob_pattern('raft_noise_correlations',
                                                 '{}_S??'.format(raft_name)),
                                    acq_jobname=acq_jobname)
    if not bias_files:
        logging.info("raft_noise_correlations: Missing bias files for raft %s",
                     raft_name)
        return None
    bias_frame_dict = defaultdict(dict)
    for item in bias_files:
        # The exposure is identified by test type, image type, and
        # seqno in the name of the folder containing the FITS files.
        seqno = os.path.dirname(item).split('_')[-1]
        # slot_name is the last '_'-delimited field before '.fits'.
        slot_name = item.split('_')[-1].split('.')[0]
        bias_frame_dict[seqno][slot_name] = item
    for seqno in bias_frame_dict:
        if len(bias_frame_dict[seqno]) == 9:
            bias_frame_files = bias_frame_dict[seqno]
            break
    logger.info("bias frame files for raft_noise_correlations:")
    for key, value in bias_frame_files.items():
        logger.info("   %s", value)
    return raft_noise_correlations(run, raft_name, bias_frame_files)


if __name__ == '__main__':
    from camera_components import camera_info
    from bot_eo_analyses import get_analysis_types, run_jh_tasks

    if 'biasnoise' in get_analysis_types():
        run_jh_tasks(read_noise_jh_task)
        run_jh_tasks(raft_jh_noise_correlations,
                     device_names=camera_info.get_raft_names())
