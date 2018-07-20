#!/usr/bin/env python
"""
Producer script for raft-level traps analysis.
"""
from eo_task_runners import run_trap_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_trap_task)
