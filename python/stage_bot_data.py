"""
Stage BOT fits files to local /scratch area for the current
harnessed job.
"""
import os
import sys
import glob
import json
import shutil
import pickle
import subprocess
import siteUtils
from camera_components import camera_info
from bot_eo_analyses import glob_pattern, bias_filename, medianed_dark_frame,\
    get_mask_files

RAFTS = set(camera_info.get_raft_names())

CCDS = set(camera_info.get_det_names())

CCD_DATA_KEYS = {'bias_frame_BOT': ('bias_frame', 'bias_stability'),
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

RAFT_DATA_KEYS = {'read_noise_BOT': ('raft_noise_correlations',),
                  'flat_pairs_BOT': ('raft_signal_correlations',),
                  'tearing_BOT': ('divisadero_tearing',)}


def write_hj_server_file(hj_fp_server_file='hj_fp_server.pkl'):
    """
    If harnessed job filepath server file is missing, query the
    eT db and create it.
    """
    if not os.path.isfile(hj_fp_server_file):
        hj_fp_server = siteUtils.HarnessedJobFilePaths()
        for analysis_type in ('badpixel', 'bias', 'dark', 'linearity',
                              'nonlinearity'):
            hj_fp_server.query_file_paths(
                siteUtils.get_analysis_run(analysis_type))
        with open(hj_fp_server_file, 'wb') as output:
            pickle.dump(hj_fp_server, output)


def get_files(data_keys, device_name):
    """
    Get the files needed by the current job for the specified device,
    e.g., 'R22_S11' or 'R22', for CCDs or rafts, respectively.
    """
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    files = set()
    for data_key in data_keys:
        pattern = glob_pattern(data_key, f'{device_name}*')
        files = files.union(
            siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname,
                                      verbose=False))
    # Remove any known bad exposures.
    bad_exposure_checker = siteUtils.BadExposureChecker()
    if not bad_exposure_checker.bad_exposures:
        return files
    good_files = set()
    for item in files:
        if bad_exposure_checker.is_bad(item):
            continue
        good_files.add(item)
    return good_files


def clean_up_scratch(run):
    """Delete scratch area associated with the specified run."""
    nodes = ('lsst-dc01 lsst-dc02 lsst-dc03 lsst-dc04 lsst-dc05 '
             'lsst-dc06 lsst-dc07 lsst-dc08 lsst-dc09 lsst-dc10').split()
    scratch_dir = os.environ.get('LCATR_SCRATCH_DIR', '/scratch')
    dest_dir = os.path.join(scratch_dir, 'bot_data', str(run))
    for node in nodes:
        cmd = ['ssh', node, 'rm', '-rf', dest_dir]
        print(' '.join(cmd))
        subprocess.check_call(cmd)


def get_isr_files(det_name, run):
    """Get bias, dark, and mask files."""
    files = set()
    try:
        bias_fn = bias_filename(run, det_name)
        if isinstance(bias_fn, str):
            files.add(bias_fn)
        else:
            files.union(bias_fn)
    except IndexError:
        pass
    try:
        files.add(medianed_dark_frame(det_name))
    except IndexError:
        pass
    files = files.union(get_mask_files(det_name))
    return files


def stage_isr_files(device_list, dest_dir):
    """
    Stage bias frame, dark frame, and mask files for the specified
    devices.
    """
    run = siteUtils.getRunNumber()
    fits_files = set()

    ccds = CCDS.intersection(device_list)
    for ccd in ccds:
        fits_files = fits_files.union(get_isr_files(ccd, run))

    rafts = RAFTS.intersection(device_list)
    for raft in rafts:
        if raft in 'R00 R04 R40 R44':
            slots = 'SG0 SG1 SW0 SW1'.split()
        else:
            slots = 'S00 S01 S02 S10 S11 S12 S20 S21 S22'.split()
        for slot in slots:
            det_name = '_'.join((raft, slot))
            fits_files = fits_files.union(get_isr_files(det_name, run))

    for src in fits_files:
        folder = os.path.basename(os.path.dirname(src))
        os.makedirs(os.path.join(dest_dir, folder), exist_ok=True)
        dest = os.path.join(dest_dir, folder, os.path.basename(src))
        if not os.path.isfile(dest):
            print('copying', src, 'to', dest)
            shutil.copy(src, dest)


def stage_files(device_list, data_keys):
    """
    Function to stage the needed raw image files from the specified
    devices (CCDs or rafts) for the current job in the scratch area.
    """
    # Gather the filenames of the needed data.
    fits_files = set()
    for device in device_list:
        fits_files = fits_files.union(get_files(data_keys, device))

    if not fits_files:
        return

    # Make the scratch directory for the BOT data.
    run_number = siteUtils.getRunNumber()
    scratch_dir = os.environ.get('LCATR_SCRATCH_DIR', '/scratch')
    dest_dir = os.path.join(scratch_dir, 'bot_data', str(run_number))
    os.makedirs(dest_dir, exist_ok=True)

    # Glob existing files to avoid re-copying or for possible clean up.
    old_files = set(glob.glob(os.path.join(dest_dir, '*', 'MC_C*.fits')))
    old_files = old_files.union(
        glob.glob(os.path.join(dest_dir, '*', 'Photodiode*.txt')))

    # Create a dict that maps src to dest file paths.  Preserve the
    # folder name of the exposure so that the PTC and flat pairs tasks
    # can identify the paired exposures.
    new_files = dict()
    frame_dirs = set()
    for src in fits_files:
        frame_dir = os.path.dirname(src)
        frame_dirs.add(frame_dir)
        folder = os.path.basename(frame_dir)
        os.makedirs(os.path.join(dest_dir, folder), exist_ok=True)
        new_files[src] = os.path.join(dest_dir, folder, os.path.basename(src))

    # Include any Photodiode_Readings*.txt files.
    for frame_dir in frame_dirs:
        for src in glob.glob(os.path.join(frame_dir, 'Photodiode*.txt')):
            new_files[src] = os.path.join(dest_dir, os.path.basename(frame_dir),
                                          os.path.basename(src))

    # Clean up unneeded files.
    unneeded_files = old_files.difference(new_files.values())
    for item in unneeded_files:
        print('removing', item)
        os.remove(item)

    # Copy the remaining files.
    for src, dest in new_files.items():
        if dest not in old_files:
            print('copying', src, 'to', dest)
            shutil.copy(src, dest)

    stage_isr_files(device_list, dest_dir)


if __name__ == '__main__':
    # Query the eT db for the archived file locations.
    write_hj_server_file()

    # Read the mapping of execution host to device lists that was
    # prepared by the TaskRunner class.
    with open('device_list_map.json', 'r') as fd:
        device_list_map = json.load(fd)
    host = sys.argv[1]
    device_list = device_list_map[host]

    job_name = os.environ['LCATR_JOB']

    ccds = CCDS.intersection(device_list)
    if ccds and job_name in CCD_DATA_KEYS:
        stage_files(ccds, CCD_DATA_KEYS[job_name])

    rafts = RAFTS.intersection(device_list)
    if rafts and job_name in RAFT_DATA_KEYS:
        stage_files(rafts, RAFT_DATA_KEYS[job_name])
