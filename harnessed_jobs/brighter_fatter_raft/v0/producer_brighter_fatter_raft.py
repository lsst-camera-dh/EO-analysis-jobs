#!/usr/bin/env python
"""
Producer script for raft-level brighter-fatter analysis.
"""

def run_bf_task(sensor_id):
    import lsst.eotest.sensor as sensorTest
    import siteUtils
    import eotestUtils

    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    flat_files = siteUtils.dependency_glob('S*/%s_flat*flat1*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('flat_pair_raft_acq'),
                                           description='Flat files:')
    bias_frame = siteUtils.dependency_glob('S*/%s_*bias*.fits' % sensor_id,
                                           jobname='fe55_raft_analysis',
                                           description='Superbias files:')[0]
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)

    task = sensorTest.BFTask()
    task.run(sensor_id, flat_files, mask_files=mask_files,
             bias_frame=bias_frame)

    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)
    siteUtils.make_png_file(plots.bf_curves,
                            '%s_brighter-fatter.png' % file_prefix,
                            bf_file='%s_bf.fits' % sensor_id)

if __name__ == '__main__':
    from multiprocessor_execution import sensor_analyses

    processes = 9                # Reserve 1 process per CCD.
    sensor_analyses(run_bf_task, processes=processes)
