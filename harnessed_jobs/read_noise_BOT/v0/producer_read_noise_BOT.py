#!/usr/bin/env python
"""
Producer script for BOT read noise analysis.
"""
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_info
from bot_eo_analyses import run_det_task_analysis, get_analysis_types, \
    raft_jh_noise_correlations

run_det_task_analysis('biasnoise')

if 'biasnoise' in get_analysis_types():
    run_device_analysis_pool(raft_jh_noise_correlations,
                             camera_info.get_raft_names(),
                             processes=None)
