#!/usr/bin/env ipython
"""
Producer script for raft-level read noise analysis.
"""
import os
from collections import defaultdict
import siteUtils
import camera_components

def run_read_noise_task(sensor_id):
    import matplotlib.pyplot as plt
    import lsst.eotest.sensor as sensorTest
    import siteUtils
    import eotestUtils
    from correlated_noise import correlated_noise

    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                           description='Fe55 bias files for read noise:')
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
    if len(bias_files) > 1:
        _, corr_fig, _ = correlated_noise(bias_files, target=0,
                                          make_plots=True, title=sensor_id)
        plt.figure(corr_fig.number)
        plt.savefig('%s_correlated_noise.png' % file_prefix)


def get_bias_files(raft_id=None):
    """Get the bias files from the Fe55 acquisition."""
    if raft_id is None:
        raft_id = os.environ['LCATR_UNIT_ID']
    raft = camera_components.Raft.create_from_etrav(raft_id)
    # Get list of bias files for each slot.
    bias_files = defaultdict(dict)
    for slot, sensor_id in raft.items():
        my_files = \
                siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                          jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                          description='Bias files for noise correlations:')
        for item in my_files:
            timestamp = item.split('_')[-1].split('.')[0]
            bias_files[timestamp][slot] = item
    for timestamp in bias_files:
        if len(bias_files[timestamp]) == 9:
            return bias_files[timestamp]
    raise RuntimeError("Could not find bias files for all nine CCDs.")


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from correlated_noise import raft_level_oscan_correlations
    from multiprocessor_execution import sensor_analyses

    processes = 9                # Reserve 1 process per CCD.
    sensor_analyses(run_read_noise_task, processes=processes)

    raft_id = os.environ['LCATR_UNIT_ID']
    run = siteUtils.getRunNumber()
    bias_files = get_bias_files(raft_id)
    title = 'Overscan correlations, {}, Run {}'.format(raft_id, run)
    plt.rcParams['figure.figsize'] = (8, 8)
    raft_level_oscan_correlations(bias_files, title=title)
    plt.savefig('{}_{}_overscan_correlations.png'.format(raft_id, run))
