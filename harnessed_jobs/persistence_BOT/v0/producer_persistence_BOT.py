#!/usr/bin/env ipython
"""
Producer script for BOT bias frame generation.
"""
import os
from persistence_jh_task import persistence_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'persistence' in get_analysis_types():
    persistence_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'persistence_BOT',
                       'v0', 'persistence_jh_task.py')
    run_python_task_or_cl_script(persistence_jh_task, persistence_task_script)
