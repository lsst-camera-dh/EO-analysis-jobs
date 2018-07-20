#!/usr/bin/env python
"""
Producer script for raft-level QE analysis.
"""
from eo_task_runners import run_qe_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_qe_task)
