#!/usr/bin/env ipython
"""
Command-line script for divisidero tearing analysis of BOT data.
"""
def raft_divisidero_tearing(raft_name):
    """JH version of divisidero tearing analysis of BOT data."""
    import json
    import matplotlib.pyplot as plt
    import siteUtils
    import lsst.eotest.raft as raftTest
    from bot_eo_analyses import get_raft_files_by_slot

    run = siteUtils.getRunNumber()

    try:
        sflat_files = get_raft_files_by_slot(raft_name, 'median_sflat.fits')
    except FileNotFoundError:
        return

    max_divisidero_tearing \
        = raftTest.ana_divisidero_tearing(sflat_files, raft_name, run)
    plt.savefig(f'{raft_name}_{run}_divisidero.png')

    with open(f'{raft_name}_{run}_max_divisidero.json', 'w') as fd:
        json.dump(max_divisidero_tearing, fd)

if __name__ == '__main__':
    import sys
    raft_name = sys.argv[1]
    raft_divisidero_tearing(raft_name)
