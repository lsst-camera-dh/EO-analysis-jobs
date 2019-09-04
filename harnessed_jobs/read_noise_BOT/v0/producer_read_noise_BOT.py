#!/usr/bin/env python
"""
Producer script for BOT read noise analysis.
"""
import os
from camera_components import camera_info
from read_noise_jh_task import read_noise_jh_task
from raft_jh_noise_correlations import raft_jh_noise_correlations
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'biasnoise' in get_analysis_types():
#    # Run the python versions of these tasks.
#    run_jh_tasks(read_noise_jh_task)
#    run_jh_tasks(raft_jh_noise_correlations,
#                 device_names=camera_info.get_raft_names())

    # Run the command-line versions.
    read_noise_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'read_noise_BOT',
                       'v0', 'read_noise_jh_task.py')
    run_jh_tasks(read_noise_jh_task_script)

    raft_jh_noise_correlations_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'read_noise_BOT',
                       'v0', 'raft_jh_noise_correlations.py')
    run_jh_tasks(raft_jh_noise_correlations_script)
