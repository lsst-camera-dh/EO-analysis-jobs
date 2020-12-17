#!/usr/bin/env ipython
"""
Script to compute signal region correlations at the raft level.
"""

def raft_jh_signal_correlations(raft_name):
    """JH version of raft-level signal-correlation analysis."""
    import os
    import logging
    import numpy as np
    import matplotlib.pyplot as plt
    import siteUtils
    from bot_eo_analyses import bias_filename, make_file_prefix, \
        append_acq_run, glob_pattern
    from signal_correlations import raft_level_signal_correlations

    logger = logging.getLogger('raft_jh_noise_correlations')
    logger.setLevel(logging.INFO)

    # Find all of the flat files, and select the pair with signal
    # levels closest to 50k e-/pixel.
    pattern = glob_pattern('flat_pairs', f'{raft_name}_*')
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    flat_files = siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname)

    target_signal = 5e4
    selected = None
    min_diff = None
    for item in flat_files:
        frame = os.path.basename(os.path.dirname(item))
        signal = frame.split('_')[3]
        diff = np.abs(target_signal - float(signal))
        if selected is None or diff < min_diff :
            selected = signal
            min_diff = diff
    flat1_files  = dict()
    flat2_files  = dict()
    for item in flat_files:
        folder = os.path.basename(os.path.dirname(item))
        basename = os.path.basename(item)
        if '_' + signal + '_' in folder:
            if raft_name in basename:
                slot = basename.split('_')[-1][:-len('.fits')]
                if 'flat0' in folder:
                    flat1_files[slot] = item
                elif 'flat1' in folder:
                    flat2_files[slot] = item
                else:
                    raise RuntimeError('unmatched flat pair file')
    logger.info('flat pair files:')
    for slot in flat1_files:
        logger.info('  ' + flat1_files[slot])
        logger.info('  ' + flat2_files[slot])
    logger.info('')

    # Find the median bias files for the target raft.
    run = siteUtils.getRunNumber()
    bias_files = dict()
    logger.info('median bias files:')
    for slot in flat1_files.keys():
        det_name = '_'.join((raft_name, slot))
        bias_files[slot] = bias_filename(run, det_name)
        logger.info('  ' + bias_files[slot])

    file_prefix = make_file_prefix(run, raft_name)
    title = append_acq_run("Imaging region correlations, "
                           f"Run {run}, {raft_name}")
    raft_level_signal_correlations(flat1_files, flat2_files, bias_files,
                                   title=title)
    plt.savefig('{}_imaging_region_correlations.png'.format(file_prefix))

if __name__ == '__main__':
    import sys
    raft_name = sys.argv[1]
    raft_jh_signal_correlations(raft_name)

