#!/usr/bin/env python
"""
Producer script for raft-level CTE analysis.
"""
from eo_task_runners import run_cte_task
from multiprocessor_execution import sensor_analyses

sensor_analyses(run_cte_task)
