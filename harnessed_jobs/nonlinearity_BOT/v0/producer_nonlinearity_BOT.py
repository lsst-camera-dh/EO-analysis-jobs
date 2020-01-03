#!/usr/bin/env ipython
"""
Producer script for BOT nonlinearity analysis.
"""
import os
from nonlinearity_jh_task import nonlinearity_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'nonlinearity' in get_analysis_types():
    nonlinearity_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                    'nonlinearity_BOT', 'v0', 'nonlinearity_jh_task.py')
    run_python_task_or_cl_script(nonlinearity_jh_task, nonlinearity_script)
