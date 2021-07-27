#!/usr/bin/env ipython
"""
Producer script for BOT dark current analysis.
"""
import os
from dark_current_jh_task import dark_current_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

analysis_types = get_analysis_types()
if 'dark' in analysis_types or 'dark_current_fit' in analysis_types:
    dark_current_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'dark_current_BOT', 'v0', 'dark_current_jh_task.py')
    run_python_task_or_cl_script(dark_current_jh_task, dark_current_script)
