#!/usr/bin/env ipython
"""
Producer script for BOT flat pairs (linearity and full-well) analysis.
"""
import os
from camera_components import camera_info
from flat_pairs_jh_task import flat_pairs_jh_task
from raft_jh_signal_correlations import raft_jh_signal_correlations
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script

if 'linearity' in get_analysis_types():
    flat_pairs_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                    'flat_pairs_BOT', 'v0', 'flat_pairs_jh_task.py')
    run_python_task_or_cl_script(flat_pairs_jh_task, flat_pairs_script)

    raft_jh_signal_correlations_script \
        =  os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                        'flat_pairs_BOT', 'v0',
                        'raft_jh_signal_correlations.py')

    # Run raft-level signal correlations just on science rafts for now.
    # Run raft-level noise correlations just on science rafts for now.
    device_names = camera_info.installed_science_rafts

    # Check if rafts are over-ridden in the lcatr.cfg file.
    override_rafts = os.environ.get('LCATR_RAFTS', None)
    if override_rafts is not None:
        device_names = [_ for _ in device_names if _[:3] in override_rafts]

    if device_names:
        run_python_task_or_cl_script(raft_jh_signal_correlations,
                                     raft_jh_signal_correlations_script,
                                     device_names=device_names)
