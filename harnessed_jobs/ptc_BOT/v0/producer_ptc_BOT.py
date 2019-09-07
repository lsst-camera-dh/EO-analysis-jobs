#!/usr/bin/env python
"""
Producer script for BOT PTC analysis.
"""
import os
from ptc_jh_task import ptc_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'ptc' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run the python version of the task.
        run_jh_tasks(ptc_jh_task)
    else:
        # Run the command-line version.
        ptc_task_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                           'harnessed_jobs', 'ptc_BOT',
                           'v0', 'ptc_jh_task.py')
        run_jh_tasks(ptc_task_script)
