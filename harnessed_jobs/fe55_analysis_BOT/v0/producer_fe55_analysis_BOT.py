#!/usr/bin/env python
"""
Producer script for BOT Fe55 analysis.
"""
import os
from fe55_jh_task import fe55_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'gain' in get_analysis_types():
    fe55_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'fe55_analysis_BOT',
                       'v0', 'fe55_jh_task.py')
    run_python_task_or_cl_script(fe55_jh_task, fe55_task_script)
