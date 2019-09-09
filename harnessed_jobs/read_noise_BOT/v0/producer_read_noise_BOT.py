#!/usr/bin/env python
"""
Producer script for BOT read noise analysis.
"""
import os
from camera_components import camera_info
from read_noise_jh_task import read_noise_jh_task
from raft_jh_noise_correlations import raft_jh_noise_correlations
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'biasnoise' in get_analysis_types():
    read_noise_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'read_noise_BOT',
                       'v0', 'read_noise_jh_task.py')
    run_python_task_or_cl_script(read_noise_jh_task, read_noise_jh_task_script)

    raft_jh_noise_correlations_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'read_noise_BOT',
                       'v0', 'raft_jh_noise_correlations.py')

    run_python_task_or_cl_script(raft_jh_noise_correlations,
                                 raft_jh_noise_correlations_script,
                                 device_names=camera_info.get_raft_names())
