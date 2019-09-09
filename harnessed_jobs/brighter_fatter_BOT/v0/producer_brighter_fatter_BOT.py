#!/usr/bin/env python
"""
Producer script for BOT brighter-fatter analysis.
"""
import os
from bf_jh_task import bf_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'brighterfatter' in get_analysis_types():
    bf_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'brighter_fatter_BOT', 'v0', 'bf_jh_task.py')
    run_python_task_or_cl_script(bf_jh_task, bf_jh_task_script)
