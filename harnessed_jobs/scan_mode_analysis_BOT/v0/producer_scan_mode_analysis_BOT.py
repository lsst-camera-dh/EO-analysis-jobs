#!/usr/bin/env python
"""
Producer script for BOT scan mode analysis.
"""
def scan_mode_analysis_jh_task(det_name):
    """JH version of scan mode analysis task."""
    import siteUtils

    run = siteUtils.getRunNumber()
    scan_mode_files = []
    return scan_mode_analysis_task(run, det_name, scan_mode_files)

from bot_eo_analyses import run_det_task_analysis
run_det_task_analysis('scan')
