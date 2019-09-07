#!/usr/bin/env python
"""
Script for BOT dark pixel analysis.
"""
def dark_defects_jh_task(det_name):
    """JH version of single sensor execution of the dark defects task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, dark_defects_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    sflat_files \
        = siteUtils.dependency_glob(glob_pattern('dark_defects', det_name),
                                    acq_jobname=acq_jobname)
    if not sflat_files:
        print("dark_defects_task: No high flux superflat files found for",
              det_name)
        return None

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))
    bias_frame = bias_filename(file_prefix)

    return dark_defects_task(run, det_name, sflat_files, mask_files=mask_files,
                             bias_frame=bias_frame)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    dark_defects_jh_task(det_name)
