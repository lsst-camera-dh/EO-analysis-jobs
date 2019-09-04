#!/usr/bin/env python
"""
Producer script for BOT tearing analysis.
"""
import os
from tearing_jh_task import tearing_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'tearing' in get_analysis_types():
#    # Run the python version.
#    run_jh_tasks(tearing_jh_task)

    # Run the command-line version.
    tearing_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'tearing_BOT', 'v0', 'tearing_jh_task.py')
    run_jh_tasks(tearing_jh_task_script)
