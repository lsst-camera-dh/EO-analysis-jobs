#!/usr/bin/env ipython
"""
Producer script for BOT flat gain stability analysis.
"""
import os
import glob
from collections import defaultdict
import siteUtils
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script
from flat_gain_stability import plot_all_rafts, plot_raft_by_amp
from flat_gain_stability_jh_task import flat_gain_stability_jh_task

if 'gainstability' in get_analysis_types():
    flat_gain_stability_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'flat_gain_stability_BOT', 'v0',
                       'flat_gain_stability_jh_task.py')
    run_python_task_or_cl_script(flat_gain_stability_jh_task,
                                 flat_gain_stability_script)

    # Make png files of plots of gains vs MJD for each raft.
    pickle_files = sorted(glob.glob('*flat_signal_sequence.pickle'))
    raft_files = defaultdict(list)
    for item in pickle_files:
        raft = os.path.basename(item)[:len('R22')]
        raft_files[raft].append(item)

    acq_run = os.environ.get('LCATR_ACQ_RUN', None)
    for raft in raft_files:
        plot_raft_by_amp(raft_files[raft], acq_run=acq_run)

    # Make focal plane summary plot of gain stability, aggregating by CCD
    plot_all_rafts(siteUtils.getRunNumber(), acq_run=acq_run)
