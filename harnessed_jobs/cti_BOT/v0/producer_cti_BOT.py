#!/usr/bin/env python
"""
Producer script for BOT PTC analysis.
"""
import os
from cte_jh_task import cte_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'cti' in get_analysis_types():
#    # Run the python version.
#    run_jh_tasks(cte_jh_task)

    # Run the command-line version.
    cte_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'cti_BOT', 'v0', 'cte_jh_task.py')
    run_jh_tasks(cte_jh_task_script)
