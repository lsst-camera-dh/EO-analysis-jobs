#!/usr/bin/env python
"""
Producer script for BOT tearing analysis.
"""
import os
from tearing_jh_task import tearing_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'tearing' in get_analysis_types():
    tearing_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'tearing_BOT', 'v0', 'tearing_jh_task.py')
    run_python_task_or_cl_script(tearing_jh_task, tearing_jh_task_script)
