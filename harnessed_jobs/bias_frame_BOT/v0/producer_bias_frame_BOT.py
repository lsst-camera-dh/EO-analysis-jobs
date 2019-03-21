#!/usr/bin/env python
"""
Producer script for BOT bias frame generation.
"""
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_info
from bot_eo_analyses import bias_frame_jh_task, get_analysis_types

# Get the list of analyses to run in this traveler from the .cfg file.
analysis_types = get_analysis_types()

# Get the list of detectors.
det_names = camera_info.get_det_names()

# Set `processes = None` to use all available cores.
processes = None

if 'bias' in analysis_types:
    run_device_analysis_pool(bias_frame_jh_task, det_names, processes=processes)
