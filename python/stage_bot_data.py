"""
Stage BOT fits files to local /scratch area for the current
harnessed job.
"""
import os
import sys
import json
import pickle
from collections import defaultdict
import siteUtils
from bot_eo_analyses import get_scan_mode_files, glob_pattern


JOB_DATA_KEYS = {'bias_frame_BOT': ('bias_frame', 'bias_stability'),
                 'fe55_analysis_BOT': ('fe55',),
                 'ptc_BOT': ('ptc',),
                 'read_noise_BOT': ('read_noise',),
                 'pixel_defects_BOT': ('bright_defects', 'dark_defects'),
                 'trap_analysis_BOT': ('traps',),
                 'persistence_BOT': ('persistence_bias', 'persistence_dark'),
                 'dark_current_BOT': ('dark_current',),
                 'flat_pairs_BOT': ('flat_pairs',),
                 'flat_gain_stability_BOT': ('tearing',),
                 'brighter_fatter_BOT': ('brighter_fatter',),
                 'cti_BOT': ('cte_high', 'cte_low'),
                 'overscan_BOT': ('overscan',),
                 'tearing_BOT': ('tearing',)}


def write_hj_server_file(hj_fp_server_file='hj_fp_server.pkl'):
    if not os.path.isfile(hj_fp_server_file):
        hj_fp_server = siteUtils.HarnessedJobFilePaths()
        for analysis_type in ('badpixel', 'bias', 'dark', 'linearity',
                              'nonlinearity'):
            hj_fp_server.query_file_paths(
                siteUtils.get_analysis_run(analysis_type))
        with open(hj_fp_server_file, 'wb') as output:
            pickle.dump(hj_fp_server, output)


def get_det_files(det_name):
    job_name = os.environ['LCATR_JOB']
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    files = set()
    for data_key in JOB_DATA_KEYS[job_name]:
        pattern = glob_pattern(data_key, det_name)
        files = files.union(
            siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname))
    return files


if __name__ == '__main__':
    # Query the eT db for the archived file locations.
    write_hj_server_file()

    # Read the mapping of execution host to device lists that was
    # prepared by the TaskRunner class.
    with open('device_map.json', 'r') as fd:
        device_map = json.load(fd)

    # Gather the filenames of the needed data.
    host = sys.argv[1]
    fits_files = set()
    for device in device_map[host]:
        fits_files.union(get_det_files(device))

    # Make the scratch directory for the BOT data.
    run_number = siteUtils.getRunNumber()
    dest_dir = f'/scratch/bot_data/{run_number}'
    os.makedirs(dest_dir, exist_ok=True)

    # Glob existing files to avoid re-copying or for possible clean up.
    old_files = set(glob.glob(os.path.join(dest_dir, 'MC_C*.fits')))

    # dict mapping src to dest file paths.
    new_files = [src: os.path.join(dest_dir, os.path.basename(src))
                 for src in fits_files]

    # Clean up unneeded files.
    unneeded_files = old_files.difference(new_files.values())
    for item in unneeded_files:
        os.remove(item)

    # Copy the remaining needed files.
    for src, dest in new_files.items():
        if dest not in old_files:
            shutil.copy(src, dest)
