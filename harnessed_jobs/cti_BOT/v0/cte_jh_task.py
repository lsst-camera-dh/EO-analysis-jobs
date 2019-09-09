#!/usr/bin/env python
"""
Script for BOT CTE analysis.
"""
def cte_jh_task(det_name):
    """JH version of single sensor execution of the CTE task."""
    import os
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, glob_pattern,\
        get_amplifier_gains, bias_filename, cte_task, plot_cte_results

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)
    acq_jobname = siteUtils.getProcessName('BOT_acq')

    sflat_high_files \
        = siteUtils.dependency_glob(glob_pattern('cte_high', det_name),
                                    acq_jobname=acq_jobname)
    sflat_low_files \
        = siteUtils.dependency_glob(glob_pattern('cte_low', det_name),
                                    acq_jobname=acq_jobname)
    if not sflat_high_files and not sflat_low_files:
        print("cte_task: Superflat files not found for detector", det_name)
        return None

    mask_files = sorted(glob.glob('{}_*mask.fits'.format(file_prefix)))

    eotest_results_file = '{}_eotest_results.fits'.format(file_prefix)
    gains = get_amplifier_gains(eotest_results_file)

    bias_frame = bias_filename(file_prefix)

    # Omit rolloff defects mask since it would mask some of the edges used
    # in the eper method.
    mask_files \
        = [item for item in mask_files if item.find('edge_rolloff') == -1]

    png_files = []
    for flux_level, sflat_files in zip(('high', 'low'),
                                       (sflat_high_files, sflat_low_files)):
        superflat_file = cte_task(run, det_name, sflat_files, gains,
                                  mask_files=mask_files, flux_level=flux_level,
                                  bias_frame=bias_frame)

        png_files.extend(plot_cte_results(run, det_name, superflat_file,
                                          eotest_results_file,
                                          mask_files=mask_files))

    png_file_list = '{}_cte_task_png_files.txt'.format(det_name)
    with open(png_file_list, 'w') as output:
        for item in png_files:
            if os.path.isfile(item):
                output.write('{}\n'.format(item))

    return None


if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    cte_jh_task(det_name)