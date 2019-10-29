#!/usr/bin/env python
"""
Producer script for BOT flat gain stability analysis.
"""
import os
from flat_gain_stability_jh_task import flat_gain_stability_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'gain_stability' in get_analysis_types():
    flat_gain_stability_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'flat_gain_stability_BOT', 'v0',
                       'flat_gain_stability_jh_task.py')
    run_python_task_or_cl_script(flat_gain_stability_jh_task,
                                 flat_gain_stability_script)
