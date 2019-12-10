#!/usr/bin/env ipython
"""
Producer script for raft-level dark defects analysis.
"""

def run_dark_pixels_task(sensor_id):
    "Single sensor execution of the dark pixels task."
    import lsst.eotest.sensor as sensorTest
    import siteUtils
    import eotestUtils

    acq_jobname = siteUtils.getProcessName('sflat_raft_acq')
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    sflat_files = siteUtils.dependency_glob('S*/%s_sflat_500_flat_H*.fits' % sensor_id,
                                            jobname=acq_jobname,
                                            description='Superflat files:')

    bias_files = siteUtils.dependency_glob('S*/%s_sflat_bias*.fits' % sensor_id,
                                           jobname=acq_jobname,
                                           description='Bias files:')
    bias_frame = eotestUtils.make_median_bias_frame(bias_files, sensor_id,
                                                    'sflat_raft_acq')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)

    task = sensorTest.DarkPixelsTask()
    task.run(sensor_id, sflat_files, mask_files, bias_frame=bias_frame)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_superflat_dark_defects.png' % file_prefix,
                            '%s_median_sflat.fits' % sensor_id,
                            title='%s, superflat for dark defects analysis' % sensor_id,
                            annotation='ADU/pixel', flatten=True, binsize=4)

if __name__ == '__main__':
    from multiprocessor_execution import sensor_analyses

    processes = 9                # Reserve 1 process per CCD.
    sensor_analyses(run_dark_pixels_task, processes=processes)
