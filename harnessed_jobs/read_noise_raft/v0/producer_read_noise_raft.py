#!/usr/bin/env python
"""
Producer script for raft-level read noise analysis.
"""
from __future__ import print_function
import multiprocessing
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
import camera_components

def run_read_noise_task(sensor_id):
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_fe55_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                           description='Fe55 files for read noise:')
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    system_noise = None

    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)

    task = sensorTest.ReadNoiseTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, bias_files, gains, system_noise=system_noise,
             mask_files=mask_files, use_overscan=True)

    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    siteUtils.make_png_file(plots.noise, '%s_noise.png' % sensor_id)

if __name__ == '__main__':
    raft_id = siteUtils.getUnitId()
    raft = camera_components.Raft.create_from_etrav(raft_id)

    # Nominally use N-1 cores in the processing pool, but ensure the
    # number is not less than 1.
    processes = max(1, multiprocessing.cpu_count() - 1)
    pool = multiprocessing.Pool(processes=processes)
    results = [pool.apply_async(run_read_noise_task, sensor_id)
               for sensor_id in raft.sensor_names]
    pool.close()
    pool.join()
    for res in results:
        res.get()
