#!/usr/bin/env python
"""
Producer script for BOT Fe55 analysis.
"""
import os
from fe55_jh_task import fe55_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'gain' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run the python version of the task.
        run_jh_tasks(fe55_jh_task)
    else:
        # Run the command-line version.
        fe55_task_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                           'harnessed_jobs', 'fe55_analysis_BOT',
                           'v0', 'fe55_jh_task.py')
        run_jh_tasks(fe55_task_script)
