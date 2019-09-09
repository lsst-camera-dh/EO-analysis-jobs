#!/usr/bin/env python
"""
Producer script for BOT bias frame generation.
"""
import os
from bias_frame_jh_task import bias_frame_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'bias' in get_analysis_types():
    bias_frame_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'bias_frame_BOT',
                       'v0', 'bias_frame_jh_task.py')
    run_python_task_or_cl_script(bias_frame_jh_task, bias_frame_task_script)
