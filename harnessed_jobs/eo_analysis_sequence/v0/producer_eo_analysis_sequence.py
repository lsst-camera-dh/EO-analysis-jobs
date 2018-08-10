#!/usr/bin/env python
import os
from eo_task_runners import *
import siteUtils
from multiprocessor_execution import sensor_analyses

task_runners = {
#    'bias': make_bias_frames,
    'gain': run_fe55_task,
    'read_noise': run_read_noise_task,
    'dark_current': run_dark_current_task,
    'bright_defects': run_bright_pixels_task,
    'dark_defects': run_dark_pixels_task,
    'ptc': run_ptc_task,
#    'brighter_fatter': brighter_fatter,
#    'gain_stability': gain_stability,
    'flat_pair': run_flat_pair_task,
#    'diffusion': fe55_psf,
    'cti': run_cte_task,
    'qe': run_qe_task,
#    'scan': scan_analysis
}

# Find the raft-level EO configuration file.
acq_config_file = os.path.join(os.environ['LCATR_CONFIG_DIR'], 'acq.cfg')
with open(acq_config_file, 'r') as acq_config:
    for line in acq_config:
        if line.startswith('rtmacqcfgfile'):
            eo_acq_config_file = line.split('=')[1].strip()

# Read in the eo analysis sequence.
with open(eo_acq_config_file, 'r') as eo_acq:
    task_names = [x.split('#')[0].split()[1] for x in eo_acq
                  if x.startswith('ANALYZE')]

# Create the mapping for the acquisition jobs to all point to the
# the multi_job_acq job.
job_map = {job: siteUtils.getProcessName('multi_job_acq') for job in
           ('fe55_raft_acq', 'dark_raft_acq', 'sflat_raft_acq',
            'ppump_raft_acq', 'flat_pair_raft_acq', 'qe_raft_acq')}

# Loop over the analysis sequence, using the sensor_analysis function
# to run the task runners via the multiprocessing module.
for task_name in task_names:
    try:
        task_runner = task_runners[task_name]
    except KeyError:
        # harnessed job implementation does not exist so skip.
        continue
    sensor_analyses(task_runner, **job_map)
