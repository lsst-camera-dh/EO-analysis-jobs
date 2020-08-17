#!/usr/bin/env ipython
"""
Producer script for BOT bias frame generation.
"""
import os
import glob
import pandas as pd
from camera_components import camera_info
import siteUtils
from bias_frame_jh_task import bias_frame_jh_task
from bot_eo_analyses import get_analysis_types, run_python_task_or_cl_script, \
    make_file_prefix

if 'bias' in get_analysis_types():
    run = siteUtils.getRunNumber()
    bias_frame_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                       'harnessed_jobs', 'bias_frame_BOT',
                       'v0', 'bias_frame_jh_task.py')
    run_python_task_or_cl_script(bias_frame_jh_task, bias_frame_task_script)

    # Combine data frames by raft.
    raft_names = camera_info.get_installed_raft_names()
    for raft_name in raft_names:
        pattern = f'{raft_name}*_bias_frame_stats.pickle'
        det_files = glob.glob(pattern)
        if not det_files:
            continue
        df = pd.concat([pd.read_pickle(_) for _ in det_files],
                       ignore_index=True)
        file_prefix = make_file_prefix(run, raft_name)
        df.to_pickle(f'{file_prefix}_bias_frame_stats.pickle')

