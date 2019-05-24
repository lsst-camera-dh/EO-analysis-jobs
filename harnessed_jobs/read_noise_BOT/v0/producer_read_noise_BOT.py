#!/usr/bin/env python
"""
Producer script for BOT read noise analysis.
"""
def read_noise_jh_task(det_name):
    """JH version of the single sensor read noise task."""
    import os
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, read_noise_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    nbias = os.environ.get('LCATR_NUM_BIAS_FRAMES', 10)

    bias_files \
        = siteUtils.dependency_glob(glob_pattern('read_noise', det_name),
                                    acq_jobname=acq_jobname)[:nbias]
    if not bias_files:
        print("read_noise_task: Needed data files are missing for detector",
              det_name)
        return None
    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    return read_noise_task(run, det_name, bias_files, gains,
                           mask_files=mask_files)


def raft_jh_noise_correlations(raft_name):
    """JH version of raft-level noise-correlation analysis."""
    import os
    import siteUtils
    from camera_components import camera_info
    from bot_eo_analyses import raft_noise_correlations, glob_pattern

    run = siteUtils.getRunNumber()
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    bias_files \
        = siteUtils.dependency_glob(glob_pattern('raft_noise_correlations',
                                                 '{}_S??'.format(raft_name)),
                                    acq_jobname=acq_jobname)
    if not bias_files:
        print("raft_noise_correlations: Missing bias files for raft",
              raft_name)
        return None
    bias_file_dict = dict()
    slot_names = camera_info.get_slot_names()
    for item in bias_files:
        for slot_name in slot_names:
            if slot_name in bias_file_dict:
                continue
            det_name = '{}_{}'.format(raft_name, slot_name)
            if det_name in os.path.basename(item):
                bias_file_dict[slot_name] = item
    return raft_noise_correlations(run, raft_name, bias_file_dict)


if __name__ == '__main__':
    from camera_components import camera_info
    from bot_eo_analyses import get_analysis_types, run_jh_tasks

    if 'biasnoise' in get_analysis_types():
        run_device_analysis_pool(read_noise_jh_task)
        run_device_analysis_pool(raft_jh_noise_correlations,
                                 device_names= camera_info.get_raft_names())
