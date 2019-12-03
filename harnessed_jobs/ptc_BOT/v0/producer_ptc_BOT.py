#!/usr/bin/env ipython
"""
Producer script for BOT PTC analysis.
"""
import os
from ptc_jh_task import ptc_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'ptc' in get_analysis_types():
    ptc_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'ptc_BOT', 'v0', 'ptc_jh_task.py')

    run_python_task_or_cl_script(ptc_jh_task, ptc_task_script)
