#!/usr/bin/env python
"""
Producer script for BOT analyses.
"""
from __future__ import print_function
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_info
from bot_eo_analyses import fe55_task, read_noise_task, dark_current_task, \
    bright_defects_task, dark_defects_task, ptc_task, flat_pairs_task, \
    cte_task, tearing_task, get_analysis_types

# Use the all of the available cores for processing.
processes = None

task_mapping = {'gain': (fe55_task,),
                'biasnoise': (read_noise_task,),
                'dark': (dark_current_task,),
                'badpixel': (bright_defects_task, dark_defects_task),
                'ptc': (ptc_task,),
                'linearity': (flat_pairs_task,),
                'cti': (cte_task,),
                'tearing': (tearing_task,)}

analysis_types = get_analysis_types()

# Detector-level analyses
det_names = camera_info.get_det_names()
for analysis_type in analysis_types:
    print("**************************************")
    print("Running analysis type %s" % analysis_type)
    print("**************************************")
    if analysis_type not in task_mapping:
        print("   not in task_mapping. skipping")
        continue
    tasks = task_mapping[analysis_type]
    for task in tasks:
        run_device_analysis_pool(task, det_names, processes=processes)

# Raft-level analyses
raft_names = camera_info.get_raft_names()
if 'biasnoise' in analysis_types:
    print("**************************************")
    print("Running analysis type biasnoise")
    print("**************************************")
    run_device_analysis_pool(raft_noise_correlations, raft_names,
                             processes=processes)

print("**************************************")
print("Running raft_results_task")
print("**************************************")
run_device_analysis_pool(raft_results_task, raft_names, processes=processes)
