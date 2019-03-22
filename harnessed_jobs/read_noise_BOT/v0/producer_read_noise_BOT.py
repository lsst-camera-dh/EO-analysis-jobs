#!/usr/bin/env python
"""
Producer script for BOT PTC analysis.
"""
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_model
from bot_eo_analyses import run_det_task_analysis, raft_jh_noise_correlations, \
    get_analysis_types

run_det_task_analysis('biasnoise')

if 'biasnoise' in get_analysis_types():
    run_device_analysis_pool(raft_jh_noise_correlations,
                             camera_model.get_raft_names(),
                             processes=None)
