#!/usr/bin/env python
"""
Producer script for raft-level dark current analysis.
"""
from eo_task_runners import run_dark_current_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_dark_current_task)
