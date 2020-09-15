"""
Example script for finding BOT data at SLAC and running some of
the BOT EO harnessed job tasks.
"""
import os
import glob
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
import bot_eo_analyses as bot_eo
from camera_components import camera_info
import siteUtils

Run = '12543'         # C-protocol run (lacks Fe55)
Fe55Run = '12535'     # 15 min Fe55 run.

raft = 'R22'
slot = 'S11'

# Harnessed job file path server.
fp_server = siteUtils.HarnessedJobFilePaths()

# Set the default run.
fp_server.acq_run = Run

# Query for the files for each run in order to cache results locally.
fp_server.query_file_paths(Run)
fp_server.query_file_paths(Fe55Run)

det_name = '_'.join((raft, slot))
# eotest results files for the different runs.
fe55_results_file = '_'.join((det_name, Fe55Run, 'eotest_results.fits'))
eotest_results_file = '_'.join((det_name, Run, 'eotest_results.fits'))

fe55_files = fp_server.get_files('BOT_acq', f'fe55_flat_*/*{det_name}.fits',
                                 run=Fe55Run)
bias_files = fp_server.get_files('BOT_acq', f'bias_bias_*/*{det_name}.fits')

# Science raft slot names
slots = [_ for _ in camera_info.get_slot_names() if _ not in 'SG0 SG1 SW0 SW1']

# Bias frames for the C-protocol run
bias_frames = sorted([os.path.basename(os.path.dirname(_)) for _ in bias_files])

# This dict is needed by the raft noise correlations task.
bias_file_dict = {_: fp_server.get_files('BOT_acq', f'{bias_frames[0]}/*{raft}_{_}.fits')[0]
                  for _ in slots}

# Find other types of files in the C-protocol run.
dark_files = fp_server.get_files('BOT_acq', f'dark_dark_*/*{det_name}.fits')
sflat_files = fp_server.get_files('BOT_acq', f'sflat_flat_*_H_*/*{det_name}.fits')
flat1_files = fp_server.get_files('BOT_acq', f'flat_*_flat1_*/*{det_name}.fits')
lambda_files = fp_server.get_files('BOT_acq', f'lambda_flat_*/*{det_name}.fits')

# Make a medianed bias frame to use in the various tasks.
bias_frame = bot_eo.make_bias_filename(Run, det_name)
bot_eo.bias_frame_task(Run, det_name, bias_files, bias_frame=bias_frame)

bot_eo.fe55_task(Fe55Run, det_name, fe55_files, bias_frame=bias_frame)
plt.close('all')    # This is needed to recover memory from matplotlib.

# Get the Fe55 gains from the fe55_results_file.  Note that this
# job generates the edge rolloff mask.
fe55_results = sensorTest.EOTestResults(fe55_results_file)
gains = dict(zip(fe55_results['AMP'], fe55_results['GAIN']))
print(gains)

# Do the read noise analysis, but only consider the first 5 bias frame
# files.  At least two are needed for this task.
bot_eo.read_noise_task(Run, det_name, bias_files[:5], gains)
plt.close('all')

bot_eo.raft_noise_correlations(Run, raft, bias_file_dict)
plt.close('all')

# Get the edge rolloff mask.
mask_files = sorted(glob.glob('_'.join((det_name, Run, '*mask.fits'))))

# The two defects tasks generate mask files, so the mask files
# need to be globbed after running the next one.
bot_eo.bright_defects_task(Run, det_name, dark_files, gains=gains,
                           mask_files=mask_files, bias_frame=bias_frame)
plt.close('all')
mask_files = sorted(glob.glob('_'.join((det_name, Run, '*mask.fits'))))

bot_eo.dark_defects_task(Run, det_name, sflat_files, mask_files=mask_files,
                         bias_frame=bias_frame)
plt.close('all')
mask_files = sorted(glob.glob('_'.join((det_name, Run, '*mask.fits'))))

dark_curr_pixels, dark95s \
    = bot_eo.dark_current_task(Run, det_name, dark_files, gains,
                               mask_files=mask_files, bias_frame=bias_frame)
bot_eo.plot_ccd_total_noise(Run, det_name, dark_curr_pixels, dark95s,
                            eotest_results_file)
plt.close('all')

superflat_file = bot_eo.cte_task(Run, det_name, sflat_files, gains,
                                 mask_files=mask_files, bias_frame=bias_frame)
bot_eo.plot_cte_results(Run, det_name, superflat_file, eotest_results_file,
                        mask_files=mask_files);
plt.close('all')

dark_frame = f'{det_name}_{Run}_median_dark_bp.fits'
bot_eo.flat_pairs_task(Run, det_name, flat1_files, gains,
                       mask_files=mask_files, bias_frame=bias_frame,
                       dark_frame=dark_frame, mondiode_func=bot_eo.mondiode_value)
plt.close('all')

bot_eo.ptc_task(Run, det_name, flat1_files, gains, mask_files=mask_files,
                bias_frame=bias_frame)
plt.close('all')

bot_eo.tearing_task(Run, det_name, flat1_files)
