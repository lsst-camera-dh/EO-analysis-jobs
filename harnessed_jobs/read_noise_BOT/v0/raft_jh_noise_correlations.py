#!/usr/bin/env python
"""
Producer script for BOT read noise analysis.
"""
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
    import sys
    raft_name = sys.argv[1]
    raft_jh_noise_correlations(raft_name)
