#!/usr/bin/env python
"""
Producer script for BOT scan mode analysis.
"""
import os
from camera_components import camera_info
from scan_mode_analysis_jh_task import scan_mode_analysis_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'scan' in get_analysis_types():
#    # Run the python version of the task.
#    run_jh_tasks(scan_mode_analysis_jh_task,
#                 device_names=camera_info.get_raft_names())

    # Run the command-line version.
    scan_mode_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'scan_mode_analysis_BOT',
                       'v0', 'scan_mode_analysis_jh_task.py')
    run_jh_tasks(scan_mode_script,
                 device_names=camera_info.get_raft_names())
