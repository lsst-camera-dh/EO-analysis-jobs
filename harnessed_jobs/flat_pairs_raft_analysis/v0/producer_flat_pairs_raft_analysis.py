#!/usr/bin/env ipython
"""
Producer script for raft-level flat pairs analysis.
"""

def run_flat_pair_task(sensor_id):
    import lsst.eotest.sensor as sensorTest
    import siteUtils
    import eotestUtils

    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    flat_files = siteUtils.dependency_glob('S*/%s_flat*flat?_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('flat_pair_raft_acq'),
                                           description='Flat files:')
    bias_frame = siteUtils.dependency_glob('%s_sflat*median_bias.fits'
                                           % sensor_id,
                                           description='Super bias frame:')[0]
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    use_exptime = True
    if siteUtils.getSiteName() == 'SLAC':
        # Since slit-width can be set individually for each exposure
        # on TS-8 at IR-2 (LSSTTD-1231), we need to use the MONDIODE
        # keyword for computing the integrated incident flux.
        use_exptime = False

    task = sensorTest.FlatPairTask()
    task.run(sensor_id, flat_files, mask_files, gains,
             linearity_spec_range=(1e4, 9e4), use_exptime=use_exptime,
             bias_frame=bias_frame)

    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    Ne_bounds = (1e4, 9e4)

    detresp_file = '%s_det_response.fits' % sensor_id
    siteUtils.make_png_file(plots.linearity,
                            '%s_linearity.png' % file_prefix,
                            detresp_file=detresp_file, max_dev=0.03,
                            use_exptime=use_exptime, Ne_bounds=Ne_bounds)
    siteUtils.make_png_file(plots.linearity_resids,
                            '%s_linearity_resids.png' % file_prefix,
                            detresp_file=detresp_file, max_dev=0.03,
                            Ne_bounds=Ne_bounds, use_exptime=use_exptime)

if __name__ == '__main__':
    from multiprocessor_execution import sensor_analyses

    processes = 9                # Reserve 1 process per CCD.
    sensor_analyses(run_flat_pair_task, processes=processes)
