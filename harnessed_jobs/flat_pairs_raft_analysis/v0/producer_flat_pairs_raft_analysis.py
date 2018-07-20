#!/usr/bin/env python
"""
Producer script for raft-level flat pairs analysis.
"""
from eo_task_runners import run_flat_pair_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_flat_pair_task)
