#!/usr/bin/env python
"""
Producer script for BOT trap analysis.
"""
import os
from traps_jh_task import traps_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'traps' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run the python version
        run_jh_tasks(traps_jh_task)
    else:
        # Run the command-line version.
        traps_jh_task_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                           'trap_analysis_BOT', 'v0', 'traps_jh_task.py')
        run_jh_tasks(traps_jh_task_script)
