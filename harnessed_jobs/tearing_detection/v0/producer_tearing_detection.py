#!/usr/bin/env python
"""
Producer script for raft-level flat pairs analysis.
"""
from __future__ import print_function
import pickle
import siteUtils
from multiprocessor_execution import sensor_analyses
from tearing_detection import tearing_detection

acq_jobs = {('flat_pair_raft_acq', 'N/A'): 'S*/%s_flat*flat?_*.fits',
            ('qe_raft_acq', 'N/A'): 'S*/%s_lambda_flat_*.fits',
            ('sflat_raft_acq', 'low_flux'): 'S*/%s_sflat_500_flat_L*.fits',
            ('sflat_raft_acq', 'high_flux'): 'S*/%s_sflat_500_flat_H*.fits'}

def run_tearing_detection(sensor_id):
    """
    Loop over the acquisition jobs and perform tearing analysis on each.
    """
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    tearing_stats = []
    for job_key, pattern in acq_jobs.items():
        job_name, subset = job_key
        flats = siteUtils.dependency_glob(pattern % sensor_id,
                                          jobname=siteUtils.getProcessName(job_name),
                                          description='Flat files:')
        tearing_found, _ = tearing_detection(flats)
        tearing_stats.append((job_name, subset, sensor_id, len(tearing_found)))

    with open('%s_tearing_stats.pkl' % file_prefix, 'wb') as output:
        pickle.dump(tearing_stats, output)

if __name__ == '__main__':
    sensor_analyses(run_tearing_detection)
