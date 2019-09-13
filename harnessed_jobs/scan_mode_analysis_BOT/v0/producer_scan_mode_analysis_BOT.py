#!/usr/bin/env python
"""
Producer script for BOT scan mode analysis.
"""
import os
from camera_components import camera_info
from scan_mode_analysis_jh_task import scan_mode_analysis_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'scan' in get_analysis_types():
    scan_mode_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'scan_mode_analysis_BOT',
                       'v0', 'scan_mode_analysis_jh_task.py')
    installed_rafts = camera_info.get_installed_raft_names()
    run_python_task_or_cl_script(scan_mode_analysis_jh_task,
                                 scan_mode_script,
                                 device_names=installed_rafts)
