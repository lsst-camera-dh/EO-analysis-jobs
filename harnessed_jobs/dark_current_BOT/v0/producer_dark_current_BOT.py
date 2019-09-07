#!/usr/bin/env python
"""
Producer script for BOT dark current analysis.
"""
import os
from dark_current_jh_task import dark_current_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'dark' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run the python version.
        run_jh_tasks(dark_current_jh_task)
    else:
        # Run the command-line version.
        dark_current_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                           'dark_current_BOT', 'v0', 'dark_current_jh_task.py')
        run_jh_tasks(dark_current_script)
