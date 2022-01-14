#!/usr/bin/env ipython
import os
import glob
import subprocess
from collections import defaultdict
from astropy.io import fits
from lsst.daf.butler import Butler
from lsst.obs.lsst import LsstCam
from bot_eo_analyses import ptc_task, make_rolloff_mask

# This needs to be set, but it's not used.
os.environ['LCATR_CONFIG_DIR'] = '.'
# Similarly, this file needs to exist, but can just be empty.
with open('acq.cfg', 'a') as fp:
    pass

CAMERA = LsstCam().getCamera()

DETECTOR = {det.getName(): detnum for detnum, det in enumerate(CAMERA)}

def get_flat_pairs(butler, run, det_name, collections=('LSSTCam/raw/all',),
                   image_type='flat', test_type='flat', staging_dir=None,
                   DETECTOR=DETECTOR):
    """
    Use the butler to find the flat pairs from a run and optionally stage
    them in the specified directory with symlinked names that conform to
    the BOT EO naming conventions so that the ptcTask.py code can be
    run on them.
    """
    raft, sensor = det_name.split('_')
    detector = DETECTOR[det_name]
    where = (f"exposure.observation_type='{image_type}' and "
             f"exposure.observation_reason='{test_type}' and "
             f"exposure.science_program='{run}' and "
             f"detector={detector}")
    dsrefs = list(butler.registry.queryDatasets('raw', instrument='LSSTCam',
                                                where=where,
                                                collections=collections))
    flat_pairs = defaultdict(list)
    for dsref in dsrefs:
        file_path = butler.getURI('raw', dsref.dataId,
                                  collections=collections).path
        if det_name not in file_path:
            continue
        with fits.open(file_path) as hdus:
            filter1 = hdus[0].header['FILTER']
            filter2 = hdus[0].header['FILTER2']
            exptime = str(hdus[0].header['EXPTIME'])
            key = '_'.join((filter1, filter2, exptime))
            flat_pairs[key].append(file_path)

    if staging_dir is None:
        return flat_pairs

    # Create symlinks in staging_dir to mimic BOT_acq `*flat[01]*` filepaths.
    os.makedirs(staging_dir, exist_ok=True)
    flat1_files = []
    for key, value in flat_pairs.items():
        if len(value) != 2:
            continue
        flat0, flat1 = value
        flat0_path = os.path.join(staging_dir,
                                  f'{det_name}_{run}_{key}_flat0.fits')

        if not os.path.islink(flat0_path):
            os.symlink(flat0, flat0_path)
        flat1_path = os.path.join(staging_dir,
                                  f'{det_name}_{run}_{key}_flat1.fits')
        if not os.path.islink(flat1_path):
            os.symlink(flat1, flat1_path)
        flat1_files.append(flat1_path)

    return flat1_files


def find_flat2(flat1):
    return flat1.replace('flat1', 'flat0')


class BiasFrame:
    def __init__(self, bias_frame_dir):
        self.bias_frame_dir = bias_frame_dir
    def __call__(self, det_name):
        if self.bias_frame_dir is None:
            return None
        def glob_file(ext):
            return glob.glob(os.path.join(self.bias_frame_dir,
                                          f'{det_name}_*_pca_bias.{ext}'))[0]
        return (glob_file('pickle'), glob_file('fits'))


if __name__ == '__main__':
    import argparse
    import multiprocessing

    parser = argparse.ArgumentParser()
    parser.add_argument('repo', type=str, help='Gen3 repo with flat pair data')
    parser.add_argument('run', type=str, help='BOT run number')
    parser.add_argument('--processes', type=int, default=1,
                        help='Number of concurrent processes. Default: 1')
    parser.add_argument('--staging_dir', type=int, default=None,
                        help=('Name of staging directory for flat pairs. '
                              'If None (default), f"{run}_flat_pairs" '
                              'will be used.'))
    parser.add_argument('--sensors', type=str, nargs='+', default=None,
                        help=('List of sensors to process, e.g., '
                              '"R22_S11 R22_S10". If None, then process all '
                              'sensors in focal plane.'))
    parser.add_argument('--bias_frame_dir', type=str, default=None,
                        help='Folder with PCA-based bias frames.')
    args = parser.parse_args()

    bias_frame = BiasFrame(args.bias_frame_dir)

    butler = Butler(args.repo)

    staging_dir = f'{args.run}_flat_pairs' if args.staging_dir is None \
                  else args.staging_dir

    gains = {}  # not actually used by ptc_task

    with multiprocessing.Pool(processes=args.processes) as pool:
        workers = []
        for det in CAMERA:
            det_name = det.getName()
            if args.sensors is not None and det_name not in args.sensors:
                continue
            flat1_files = get_flat_pairs(butler, args.run, det_name,
                                         staging_dir=staging_dir)
            mask_files = (make_rolloff_mask(args.run, det_name, flat1_files[0],
                                            outdir=staging_dir),)
            pars = args.run, det_name, flat1_files, gains
            kwds = dict(flat2_finder=find_flat2, mask_files=mask_files,
                        bias_frame=bias_frame(det_name))
            workers.append(pool.apply_async(ptc_task, pars, kwds))
        pool.close()
        pool.join()
        _ = [worker.get() for worker in workers]

    # Clean up staging_dir
    subprocess.check_call(f'rm -rf {staging_dir}', shell=True)
