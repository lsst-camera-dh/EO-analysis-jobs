#!/usr/bin/env python
"""
Producer script for BOT tearing analysis.
"""
import os
from camera_components import camera_info
from tearing_jh_task import tearing_jh_task
from raft_divisidero_tearing import raft_divisidero_tearing
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'tearing' in get_analysis_types():
    tearing_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'tearing_BOT', 'v0', 'tearing_jh_task.py')
    run_python_task_or_cl_script(tearing_jh_task, tearing_jh_task_script)

    divisidero_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'tearing_BOT', 'v0', 'raft_divisidero_tearing.py')
    installed_rafts = camera_info.get_installed_raft_names()
    run_python_task_or_cl_script(raft_divisidero_tearing,
                                 divisidero_jh_task_script,
                                 device_names=installed_rafts)
