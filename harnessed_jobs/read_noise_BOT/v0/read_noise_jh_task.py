#!/usr/bin/env ipython
"""
Producer script for BOT read noise analysis.
"""
def read_noise_jh_task(det_name):
    """JH version of the single sensor read noise task."""
    import os
    import glob
    import logging
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, read_noise_task, get_mask_files

    logger = logging.getLogger('read_noise_jh_task')
    logger.setLevel(logging.INFO)

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    nbias = os.environ.get('LCATR_NUM_BIAS_FRAMES', 10)

    bias_files \
        = siteUtils.dependency_glob(glob_pattern('read_noise', det_name),
                                    acq_jobname=acq_jobname)[:nbias]
    if not bias_files:
        logger.info("read_noise_task: Needed data files are missing "
                    "for detector %s", det_name)
        return None
    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    mask_files = get_mask_files(det_name)

    return read_noise_task(run, det_name, bias_files, gains,
                           mask_files=mask_files)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    read_noise_jh_task(det_name)
