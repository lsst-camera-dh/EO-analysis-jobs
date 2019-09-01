#!/usr/bin/env python
"""
Producer script for BOT bias frame generation.
"""
from bias_frame_jh_task import bias_frame_jh_task

if __name__ == '__main__':
    import os
    import sys
    from bot_eo_analyses import get_analysis_types, run_jh_tasks
    if 'bias' in get_analysis_types():
        # Run the python version of the task.
        run_jh_tasks(bias_frame_jh_task)

#        # Run the command-line version.
#        bias_frame_task_script \
#            = os.path.join(os.path.environ['EOANALYSISJOBSDIR'],
#                           'bias_frame_BOT', 'bias_Frame_jh_task.py'))
#        run_jh_tasks(bias_frame_task_script)
