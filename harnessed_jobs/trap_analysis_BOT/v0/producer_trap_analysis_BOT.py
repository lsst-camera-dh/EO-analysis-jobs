#!/usr/bin/env python
"""
Producer script for BOT trap analysis.
"""
import os
from traps_jh_task import traps_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'traps' in get_analysis_types():
    traps_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'trap_analysis_BOT', 'v0', 'traps_jh_task.py')
    run_python_task_or_cl_script(traps_jh_task, traps_jh_task_script)
