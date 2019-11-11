#!/usr/bin/env python
"""
Producer script for BOT Fe55 analysis.
"""
import os
from fe55_jh_task import fe55_jh_task
from gain_stability_jh_task import gain_stability_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

job_dir = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'fe55_analysis_BOT', 'v0')

analysis_types = get_analysis_types()

if 'gain' in analysis_types or 'gainstability' in analysis_types:
    fe55_task_script = os.path.join(job_dir, 'fe55_jh_task.py')
    run_python_task_or_cl_script(fe55_jh_task, fe55_task_script)

if 'gainstability' in analysis_types:
    gain_stability_script = os.path.join(job_dir, 'gain_stability_jh_task.py')
    run_python_task_or_cl_script(gain_stability_jh_task, gain_stability_script)
