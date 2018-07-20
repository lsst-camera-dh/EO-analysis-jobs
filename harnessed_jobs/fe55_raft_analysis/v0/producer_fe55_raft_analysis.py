#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
from eo_task_runners import run_fe55_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_fe55_task)
