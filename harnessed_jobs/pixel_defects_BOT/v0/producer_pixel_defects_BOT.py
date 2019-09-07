#!/usr/bin/env python
"""
Producer script for BOT bright and dark pixel analyses.
"""
import os
from bright_defects_jh_task import bright_defects_jh_task
from dark_defects_jh_task import dark_defects_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'badpixel' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run python versions.
        run_jh_tasks(bright_defects_jh_task, dark_defects_jh_task)
    else:
        # Run command-line versions.
        bright_defects_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                           'harnessed_jobs', 'pixel_defects_BOT',
                           'v0', 'bright_defects_jh_task.py')
        run_jh_tasks(bright_defects_script)

        dark_defects_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                           'harnessed_jobs', 'pixel_defects_BOT',
                           'v0', 'dark_defects_jh_task.py')
        run_jh_tasks(dark_defects_script)
