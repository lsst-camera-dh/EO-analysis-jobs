#!/usr/bin/env python
"""
Producer script for raft-level read noise analysis.
"""
from eo_task_runners import run_read_noise_task, raft_overscan_correlations
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_read_noise_task)
raft_overscan_correlations()
