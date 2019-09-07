#!/usr/bin/env python
"""
Producer script for BOT flat pairs (linearity and full-well) analysis.
"""
import os
from flat_pairs_jh_task import flat_pairs_jh_task
from bot_eo_analyses import get_analysis_types, run_jh_tasks

if 'linearity' in get_analysis_types():
    if os.environ.get('LCATR_USE_PARSL', False) != 'True':
        # Run the python version.
        run_jh_tasks(flat_pairs_jh_task)
    else:
        # Run the command-line version.
        flat_pairs_script \
            = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                           'flat_pairs_BOT', 'v0', 'flat_pairs_jh_task.py')
        run_jh_tasks(flat_pairs_script)
