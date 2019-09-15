#!/usr/bin/env python
"""
Producer script for raft-level Fe55 analysis.
"""
# Note that parsl python_apps require imports to be within the bodies
# of these functions since it is the serialized versions of these
# objects that are executed remotely.

def run_fe55_task(sensor_id):
    "Single sensor execution of the Fe55 analysis task."
    import os
    import glob
    import lsst.eotest.image_utils as imutils
    import lsst.eotest.sensor as sensorTest
    import siteUtils
    import eotestUtils

    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    acq_jobname = siteUtils.getProcessName('fe55_raft_acq')
    fe55_files = siteUtils.dependency_glob('S*/%s_fe55_fe55_*.fits' % sensor_id,
                                           jobname=acq_jobname,
                                           description='Fe55 files:')
    # Reverse sort the fe55 files to avoid transient effects arising
    # from using frames taken right after cool down that could bias
    # the gain measurement.
    fe55_files = sorted(fe55_files, reverse=True)
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                           jobname=acq_jobname,
                                           description='Bias files:')
    bias_frame = eotestUtils.make_median_bias_frame(bias_files, sensor_id,
                                                    'fe55_raft_acq', skip=0)

    #
    # Create a png zoom of the upper right corner of segment 1 for an Fe55
    # exposure for inclusion in the test report.
    #
    print("processing fe55_zoom:", fe55_files[0])
    siteUtils.make_png_file(sensorTest.fe55_zoom,
                            '%(file_prefix)s_fe55_zoom.png' % locals(),
                            fe55_files[0], size=250, amp=1,
                            annotation='ADU/pixel')

    #
    # Perform analysis of 9-pixel statistics for Fe55 charge clusters.
    #
    try:
        pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=sensor_id)

        siteUtils.make_png_file(pixel_stats.pixel_hists,
                                '%s_fe55_p3_p5_hists.png' % file_prefix,
                                pix0='p3', pix1='p5')

        siteUtils.make_png_file(pixel_stats.pixel_diff_profile,
                                '%s_fe55_p3_p5_profiles.png' % file_prefix,
                                pixel_coord='x', pix0='p3', pix1='p5')

        siteUtils.make_png_file(pixel_stats.apflux_profile,
                                '%s_fe55_apflux_serial.png' % file_prefix)

        siteUtils.make_png_file(pixel_stats.apflux_profile,
                                '%s_fe55_apflux_parallel.png' % file_prefix,
                                pixel_coord='y')

    except Exception as eobj:
        print("Exception raised while creating pixel statistics plots:")
        print(str(eobj))
        print("Skipping these plots.")

    # Roll-off defects mask needs an input file to get the vendor
    # geometry, and will be used for all analyses.
    rolloff_mask_file = '%s_rolloff_defects_mask.fits' % sensor_id
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    task = sensorTest.Fe55Task()
    task.config.temp_set_point = -100.
    hist_nsig = 10
    task.run(sensor_id, fe55_files, (rolloff_mask_file,), bias_frame=bias_frame,
             accuracy_req=0.01, hist_nsig=hist_nsig)

    # Fe55 gain and psf analysis results plots for the test report.
    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    siteUtils.make_png_file(plots.gains,
                            '%s_gains.png' % file_prefix)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_median_bias.png' % file_prefix,
                            bias_frame,
                            title='%s, median bias frame' % sensor_id,
                            annotation='ADU/pixel, overscan-subtracted')

    fe55_file = glob.glob('%s_psf_results*.fits' % sensor_id)[0]
    siteUtils.make_png_file(plots.fe55_dists,
                            '%s_fe55_dists.png' % file_prefix,
                            fe55_file=fe55_file, xrange_scale=3)

    siteUtils.make_png_file(plots.psf_dists,
                            '%s_psf_dists.png' % file_prefix,
                            fe55_file=fe55_file)

def write_nominal_gains(sensor_id, gain=1):
    """Write nominal gains in lieu of doing the Fe55 analysis."""
    import lsst.eotest.sensor as sensorTest

    results_file = '%s_eotest_results.fits' % sensor_id
    results = sensorTest.EOTestResults(results_file)
    for amp in range(1, 17):
        results.add_seg_result(amp, 'GAIN', gain)
    results.write(clobber=True)


if __name__ == '__main__':
    import os
    from multiprocessor_execution import sensor_analyses

    processes = 9     # Reserve 1 process per CCD.
    if os.environ.get("LCATR_SKIP_FE55_ANALYSIS", "False") == "True":
        sensor_analyses(write_nominal_gains, processes=processes)
    else:
        sensor_analyses(run_fe55_task, processes=processes)
