#!/usr/bin/env python
"""
Producer script for raft-level bright defects analysis.
"""
from eo_task_runners import run_bright_pixels_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_bright_pixels_task)
