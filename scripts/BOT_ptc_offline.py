import os
from collections import defaultdict
from astropy.io import fits
from lsst.daf.butler import Butler
from lsst.obs.lsst import LsstCam
from bot_eo_analyses import ptc_task


def get_flat_pairs(butler, run, det_name, collections=('LSSTCam/raw/all',),
                   image_type='flat', test_type='flat', staging_dir=None):
    raft, sensor = det_name.split('_')
    where = (f'exposure.observation_type="{image_type}" and '
             f'exposure.observation_reason="{test_type}" and '
             f'exposure.science_program in ("{run}") and '
             f'exposure.detector_group="{raft}" and '
             f'exposure.detector_name="{sensor}"')
    dsrefs = butler.registry.queryDatasets('raw', instrument='LSSTCam',
                                           where=where,
                                           collections=collections)
    flat_pairs = defaultdict(list)
    for dsref in dsrefs:
        file_path = butler.getURI(dsref.dataId, collections=collections).path
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
    flat0_files = []
    for key, value in flat_pairs.items():
        if len(value) != 2:
            continue
        flat0, flat1 = value
        flat0_path = os.path.join(staging_dir,
                                  f'{det_name}_{run}_{key}_flat0.fits')

        os.symlink(flat0, flat0_path)
        flat1_path = os.path.join(staging_dir,
                                  f'{det_name}_{run}_{key}_flat1.fits')
        os.symlink(flat1, flat1_path)
        flat0_files.append(flat0_path)

    return flat0_files


def find_flat1(flat0):
    return flat0.replace('flat0', 'flat1')


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
    args = parser.parse_args()

    camera = LsstCam().getCamera()

    butler = Butler(args.repo)

    staging_dir = f'{args.run}_flat_pairs' if args.staging_dir is None \
                  else args.staging_dir

    gains = {}  # not actually used by ptc_task

    with multiprocessing.Pool(processes=args.processes) as pool:
        workers = []
        for det in camera:
            det_name = det.getName()
            flat0_files = get_flat_pairs(butler, args.run, det_name,
                                         staging_dir=staging_dir)
            pars = args.run, det_name, flat0_files, gains
            kwds = dict(flat2_finder=find_flat1)
            workers.append(pool.apply_async(ptc_task, pars, kwds))
        pool.close()
        pool.join()
        _ = [worker.get() for worker in workers]
