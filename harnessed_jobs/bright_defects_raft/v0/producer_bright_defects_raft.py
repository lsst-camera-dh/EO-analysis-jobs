#!/usr/bin/env python
"""
Producer script for raft-level bright defects analysis.
"""

def run_bright_pixels_task(sensor_id):
    "Single sensor execution of the bright pixels task."
    import lsst.eotest.sensor as sensorTest
    import lsst.eotest.image_utils as imutils
    import siteUtils
    import eotestUtils

    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    acq_jobname = siteUtils.getProcessName('dark_raft_acq')
    dark_files = siteUtils.dependency_glob('S*/%s_dark_dark_*.fits' % sensor_id,
                                           jobname=acq_jobname,
                                           description='Dark files:')
    bias_files = siteUtils.dependency_glob('S*/%s_dark_bias_*.fits' % sensor_id,
                                           jobname=acq_jobname,
                                           description='Bias files:')
    bias_files = sorted(bias_files)[1:] # Skip the first frame
    nframes = len(bias_files)
    bias_frame = f'{sensor_id}_dark_acq_median_bias_{nframes}.fits'
    amp_geom = sensorTest.makeAmplifierGeometry(bias_files[0])
    imutils.superbias_file(bias_files, amp_geom.serial_overscan, bias_frame)

    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.BrightPixelsTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, dark_files, mask_files, gains, bias_frame=bias_frame)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_medianed_dark.png' % file_prefix,
                            '%s_median_dark_bp.fits' % sensor_id,
                            title='%s, medianed dark for bright defects analysis' % sensor_id,
                            annotation='e-/pixel, gain-corrected, bias-subtracted',
                            bias_frame=bias_frame, gains=gains, binsize=4)


if __name__ == '__main__':
    from multiprocessor_execution import sensor_analyses

    processes = 9                # Reserve 1 process per CCD.
    sensor_analyses(run_bright_pixels_task, processes=processes)
