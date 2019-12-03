#!/usr/bin/env ipython
"""
Producer script for BOT bright and dark pixel analyses.
"""
import os
from bright_defects_jh_task import bright_defects_jh_task
from dark_defects_jh_task import dark_defects_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'badpixel' in get_analysis_types():
    bright_defects_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'pixel_defects_BOT',
                       'v0', 'bright_defects_jh_task.py')
    run_python_task_or_cl_script(bright_defects_jh_task, bright_defects_script)

    dark_defects_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'pixel_defects_BOT',
                       'v0', 'dark_defects_jh_task.py')
    run_python_task_or_cl_script(dark_defects_jh_task, dark_defects_script)
