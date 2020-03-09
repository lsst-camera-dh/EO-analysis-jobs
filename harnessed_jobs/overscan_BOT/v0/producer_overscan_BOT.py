#!/usr/bin/env ipython
"""
Producer script for BOT overscan analysis.
"""
import os
from overscan_jh_task import overscan_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'overscan' in get_analysis_types():
    overscan_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'overscan_BOT', 'v0', 'overscan_jh_task.py')

    run_python_task_or_cl_script(overscan_jh_task, overscan_task_script)
