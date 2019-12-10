#!/usr/bin/env ipython
"""
Producer script for BOT PTC analysis.
"""
import os
from cte_jh_task import cte_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'cti' in get_analysis_types():
    cte_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'cti_BOT', 'v0', 'cte_jh_task.py')
    run_python_task_or_cl_script(cte_jh_task, cte_jh_task_script)
