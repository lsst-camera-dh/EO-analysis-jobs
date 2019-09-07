#!/usr/bin/env python
"""
Producer script for BOT bias frame generation.
"""
import os
from bias_frame_jh_task import bias_frame_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'bias' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run the python version of the task in a
        # multiprocessing.Pool
        run_jh_tasks(bias_frame_jh_task)
    else:
        # Run the command-line version with parsl.
        bias_frame_task_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                           'harnessed_jobs', 'bias_frame_BOT',
                           'v0', 'bias_frame_jh_task.py')
        run_jh_tasks(bias_frame_task_script)
