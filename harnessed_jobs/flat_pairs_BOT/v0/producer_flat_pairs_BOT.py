#!/usr/bin/env python
"""
Producer script for BOT flat pairs (linearity and full-well) analysis.
"""
import os
from flat_pairs_jh_task import flat_pairs_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'linearity' in get_analysis_types():
    flat_pairs_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                    'flat_pairs_BOT', 'v0', 'flat_pairs_jh_task.py')
    run_python_task_or_cl_script(flat_pairs_jh_task, flat_pairs_script)
