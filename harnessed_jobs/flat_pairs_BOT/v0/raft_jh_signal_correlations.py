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

    # Use the high signal superflat files.
    pattern = glob_pattern('raft_signal_correlations', f'{raft_name}_S??')
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    sflat_files = siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname)
    folders = sorted(list(set([os.path.basename(os.path.dirname(_))
                               for _ in sflat_files])))
    logger.info(f'folders: {folders}')

    flat1_files = dict()
    flat2_files = dict()
    for item in sflat_files:
        folder = os.path.basename(os.path.dirname(item))
        if folder not in folders[:2]:
            continue
        logger.info(f'item: {item}')
        logger.info(f'folder: {folder}')
        basename = os.path.basename(item)
        logger.info(f'basename: {basename}')
        slot = basename.split('_')[-1][:-len('.fits')]
        if folder == folders[0]:
            flat1_files[slot] = item
        elif folder == folders[1]:
            flat2_files[slot] = item
    logger.info('flat pair files:')
    for slot in flat1_files:
        logger.info('  ' + flat1_files[slot])
        logger.info('  ' + flat2_files[slot])

    # Find the median bias files for the target raft.
    run = siteUtils.getRunNumber()
    bias_files = dict()
    for slot in flat1_files.keys():
        det_name = '_'.join((raft_name, slot))
        bias_files[slot] = bias_filename(run, det_name)

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

