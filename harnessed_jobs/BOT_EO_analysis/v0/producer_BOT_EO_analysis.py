#!/usr/bin/env python
"""
Producer script for BOT analyses.
"""
from __future__ import print_function
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_info
from bot_eo_analyses import fe55_jh_task, bias_frame_jh_task, \
    read_noise_jh_task, \
    dark_current_jh_task, bright_defects_jh_task, dark_defects_jh_task, \
    ptc_jh_task, flat_pairs_jh_task, bf_jh_task, \
    cte_jh_task, tearing_jh_task, \
    get_analysis_types, raft_jh_noise_correlations, raft_results_task

# Use the all of the available cores for processing.
processes = None

task_mapping = {'gain': (fe55_jh_task,),
                'bias': (bias_frame_jh_task,),
                'biasnoise': (read_noise_jh_task,),
                'dark': (dark_current_jh_task,),
                'badpixel': (bright_defects_jh_task, dark_defects_jh_task),
                'ptc': (ptc_jh_task,),
                'brighterfatter': (bf_jh_task,),
                'linearity': (flat_pairs_jh_task,),
                'cti': (cte_jh_task,),
                'tearing': (tearing_jh_task,)}

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
    run_device_analysis_pool(raft_jh_noise_correlations, raft_names,
                             processes=processes)

print("**************************************")
print("Running raft_results_task")
print("**************************************")
run_device_analysis_pool(raft_results_task, raft_names, processes=processes)
