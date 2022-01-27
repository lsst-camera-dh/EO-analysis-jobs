#!/usr/bin/env ipython
"""
Command-line script for divisidero tearing analysis of BOT data.
"""
def raft_divisidero_tearing(raft_name):
    """JH version of divisidero tearing analysis of BOT data."""
    import os
    from collections import defaultdict
    import json
    import matplotlib.pyplot as plt
    import siteUtils
    from lsst.eotest.sensor.cteTask import superflat
    import lsst.eotest.raft as raftTest
    from bot_eo_analyses import glob_pattern, bias_filename, get_mask_files

    run = siteUtils.getRunNumber()

    pattern = glob_pattern('divisadero_tearing', f'{raft_name}_*')
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    bot_data = siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname)
    if not bot_data:
        return

    sflat_files = defaultdict(list)
    for item in bot_data:
        slot = item.split('.')[0].split('_')[-1]
        sflat_files[slot].append(item)

    median_sflats = dict()
    for slot, files in sflat_files.items():
        det_name = '_'.join((raft_name, slot))
        outfile = f'{det_name}_{run}_median_sflat.fits'
        bias_frame = bias_filename(run, det_name)
        median_sflats[slot] = superflat(files, outfile=outfile,
                                        bias_frame=bias_frame)

    mask_files = dict()
    for slot in sflat_files:
        det_name = '_'.join((raft_name, slot))
        mask_files[slot] = get_mask_files(det_name)

    title = f'Run {run} {raft_name}'
    acq_run = os.environ.get('LCATR_ACQ_RUN', None)
    if acq_run is not None:
        title += f' (acq {acq_run})'

    max_divisidero_tearing \
        = raftTest.ana_divisidero_tearing(median_sflats, mask_files,
                                          title=title)
    plt.savefig(f'{raft_name}_{run}_divisidero.png')

    with open(f'{raft_name}_{run}_max_divisidero.json', 'w') as fd:
        json.dump(max_divisidero_tearing, fd)

if __name__ == '__main__':
    import sys
    raft_name = sys.argv[1]
    raft_divisidero_tearing(raft_name)
