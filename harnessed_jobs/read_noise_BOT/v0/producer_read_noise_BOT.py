#!/usr/bin/env ipython
"""
Producer script for BOT read noise analysis.
"""
import os
from camera_components import camera_info
from read_noise_jh_task import read_noise_jh_task
from raft_jh_noise_correlations import raft_jh_noise_correlations
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'biasnoise' in get_analysis_types():
    read_noise_jh_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'read_noise_BOT',
                       'v0', 'read_noise_jh_task.py')
    run_python_task_or_cl_script(read_noise_jh_task, read_noise_jh_task_script)

    raft_jh_noise_correlations_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'read_noise_BOT',
                       'v0', 'raft_jh_noise_correlations.py')

    # Run raft-level noise correlations just on science rafts for now.
    device_names = camera_info.installed_science_rafts

    # Check if rafts are over-ridden in the lcatr.cfg file.
    override_rafts = os.environ.get('LCATR_RAFTS', None)
    if override_rafts is not None:
        device_names = [_ for _ in device_names if _[:3] in override_rafts]

    if device_names:
        run_python_task_or_cl_script(raft_jh_noise_correlations,
                                     raft_jh_noise_correlations_script,
                                     device_names=device_names)
