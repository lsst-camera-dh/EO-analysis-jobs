#!/usr/bin/env python
"""
Producer script for BOT scan mode analysis.
"""
def scan_mode_analysis_jh_task(raft_name):
    """JH version of scan mode analysis task."""
    import siteUtils
    from bot_eo_analyses import get_scan_mode_files, scan_mode_analysis_task

    run = siteUtils.getRunNumber()
    scan_mode_files = get_scan_mode_files(raft_name)
    return scan_mode_analysis_task(run, raft_name, scan_mode_files)


if __name__ == '__main__':
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'scan' in get_analysis_types():
        run_jh_tasks(scan_mode_analysis_jh_task)
