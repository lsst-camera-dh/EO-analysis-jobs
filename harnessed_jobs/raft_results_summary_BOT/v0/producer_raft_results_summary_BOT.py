#!/usr/bin/env python
"""
Producer script for BOT raft-level results summaries.
"""
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_info
from bot_eo_analyses import raft_results_task

run_device_analysis_pool(raft_results_task, camera_info.get_raft_names(),
                         processes=None)
