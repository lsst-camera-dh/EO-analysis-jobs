#!/usr/bin/env ipython
"""
Producer script for BOT Fe55 analysis.
"""
import os
import glob
from collections import defaultdict
from fe55_jh_task import fe55_jh_task
from gain_stability_jh_task import gain_stability_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script
from fe55_gain_stability import plot_raft_fe55_gains_by_amp, \
    plot_all_raft_fe55_gains

job_dir = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'fe55_analysis_BOT', 'v0')

analysis_types = get_analysis_types()

if 'gain' in analysis_types or 'gainstability' in analysis_types:
    fe55_task_script = os.path.join(job_dir, 'fe55_jh_task.py')
    run_python_task_or_cl_script(fe55_jh_task, fe55_task_script)

if 'gainstability' in analysis_types:
    gain_stability_script = os.path.join(job_dir, 'gain_stability_jh_task.py')
    run_python_task_or_cl_script(gain_stability_jh_task, gain_stability_script)

    # Make png files of plots of gains vs MJD for each raft.
    pickle_files = sorted(glob.glob('*_gain_sequence.pickle'))
    raft_files = defaultdict(list)
    for item in pickle_files:
        raft = os.path.basename(item)[:len('R22')]
        raft_files[raft].append(item)

    for raft in raft_files:
        plot_raft_fe55_gains_by_amp(raft_files[raft])

    # Make focal plane summary plot of gain stability plotting by CCD.
    plot_all_raft_fe55_gains(raft_files)
