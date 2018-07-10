#!/usr/bin/env python
"""
Producer script for raft-level read noise analysis.
"""
from __future__ import print_function
import os
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from correlated_noise import correlated_noise, raft_level_oscan_correlations
from multiprocessor_execution import sensor_analyses
import camera_components


def run_read_noise_task(sensor_id):
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
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

    # Compute amp-amp correlated noise.
    _, corr_fig, _ = correlated_noise(bias_files, target=0,
                                      make_plots=True, title=sensor_id)
    plt.figure(corr_fig.number)
    plt.savefig('%s_correlated_noise.png' % file_prefix)


def get_bias_files(raft_id=None):
    """Get the bias files from the Fe55 acquisition."""
    if raft_id is None:
        raft_id = os.environ['LCATR_UNIT_ID']
    raft = camera_components.Raft.create_from_etrav(raft_id)
    bias_files = dict()
    for slot, sensor_id in raft.items():
        bias_files[slot] \
            = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                        jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                        description='Bias files for noise correlations:')[0]
    return bias_files

if __name__ == '__main__':
    sensor_analyses(run_read_noise_task)

    raft_id = os.environ['LCATR_UNIT_ID']
    run = siteUtils.getRunNumber()
    bias_files = get_bias_files(raft_id)
    title = 'Overscan correlations, {}, Run {}'.format(raft_id, run)
    plt.rcParams['figure.figsize'] = (8, 8)
    raft_level_oscan_correlations(bias_files, title=title)
    plt.savefig('{}_{}_overscan_correlations.png'.format(raft_id, run))
