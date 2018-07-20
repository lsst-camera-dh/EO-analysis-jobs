#!/usr/bin/env python
"""
Producer script for raft-level PTC analysis.
"""
from eo_task_runner import run_ptc_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_ptc_task)
