import os
import glob
import shutil
import multiprocessing
from astropy.io import fits
import lsst.eotest.sensor as sensorTest
from bot_eo_analyses import bias_frame_task, flat_gain_stability_task, \
    make_file_prefix
from flat_gain_stability import plot_raft_by_amp

_DET_NAMES = tuple(f'R22_S{slot}' for slot in
                   '00 01 02 10 11 12 20 21 22'.split())

def find_frames(data_dir, imgtype='BIAS', testtype='FLAT',
                glob_pattern='TS_C_*/*R22_S11.fits'):
    files = sorted(glob.glob(os.path.join(data_dir, glob_pattern)))
    frame_dirs = []
    for item in files:
        with fits.open(item) as hdus:
            header = hdus[0].header
            if header['IMGTYPE'] == imgtype and header['TESTTYPE'] == testtype:
                frame_dirs.append(os.path.dirname(item))
    return frame_dirs


def make_mask_file(run, det_name, template_file):
    file_prefix = make_file_prefix(run, det_name)
    rolloff_mask_file = f'{file_prefix}_edge_rolloff_mask.fits'
    sensorTest.rolloff_mask(template_file, rolloff_mask_file)
    return rolloff_mask_file


def run_bias_frame_task(run, frame_dirs, det_names=_DET_NAMES, processes=9,
                        outdir=None, nskip=5):
    with multiprocessing.Pool(processes=processes) as pool:
        futures = []
        for det_name in det_names:
            bias_files = []
            for frame_dir in frame_dirs[nskip:]:
                frame = os.path.basename(frame_dir)
                bias_files.append(
                    os.path.join(frame_dir, f'{frame}_{det_name}.fits'))
            assert all(os.path.isfile(_) for _ in bias_files)
            args = run, det_name, bias_files
            futures.append(pool.apply_async(bias_frame_task, args))
        pool.close()
        pool.join()

        for future in futures:
            bias_frame, pca_file = future.get()

    # Move data products to output directory.
    if outdir is None:
        outdir = 'outputs/bias_models'
    os.makedirs(outdir, exist_ok=True)
    for pattern in ('*median_bias.fits', '*pca_bias.*', '*rolloff_mask.fits'):
        for src in glob.glob(pattern):
            dest = os.path.join(outdir, os.path.basename(src))
            shutil.move(src, dest)

    return outdir


def run_flat_gain_stability_task(run, frame_dirs, bias_dir,
                                 det_names=_DET_NAMES, processes=9,
                                 outdir=None):
    def add_bias_dir(x):
        return os.path.join(bias_dir, x)

    with multiprocessing.Pool(processes=processes) as pool:
        futures = []
        for det_name in det_names:
            prefix = make_file_prefix(run, det_name)
            flat_files = []
            for frame_dir in frame_dirs:
                frame = os.path.basename(frame_dir)
                flat_files.append(
                    os.path.join(frame_dir, f'{frame}_{det_name}.fits'))
            assert all(os.path.isfile(_) for _ in flat_files)
            if bias_dir is None:
                mask_files = (make_mask_file(run, det_name, flat_files[0]),)
                bias_frame = None
            else:
                mask_files = (add_bias_dir(f'{prefix}_edge_rolloff_mask.fits'),)
                bias_frame = (add_bias_dir(f'{prefix}_pca_bias.pickle'),
                              add_bias_dir(f'{prefix}_pca_bias.fits'))

            args = run, det_name, flat_files
            kwds = dict(mask_files=mask_files, bias_frame=bias_frame)
            futures.append(
                pool.apply_async(flat_gain_stability_task, args, kwds))
        pool.close()
        pool.join()

        for future in futures:
            df, outfile = future.get()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('ts8_data_dir', type=str,
                        help='directory with TS8 data')
    parser.add_argument('--processes', type=int, default=9,
                        help='number of parallel processses')
    parser.add_argument('--outdir', type=str, default='.',
                        help='output directory for results')
    args = parser.parse_args()

    Run = '-'.join(('TS8', os.path.basename(args.ts8_data_dir.rstrip('/'))))
    bias_dir = None
    bias_frame_dirs = find_frames(args.ts8_data_dir)
    if bias_frame_dirs:
        bias_dir = run_bias_frame_task(Run, bias_frame_dirs,
                                       processes=args.processes)

    flat_frame_dirs = find_frames(args.ts8_data_dir, imgtype='FLAT')
    run_flat_gain_stability_task(Run, flat_frame_dirs, bias_dir,
                                 processes=args.processes)

    raft_files = sorted(glob.glob('*flat_signal_sequence.pickle'))
    plot_raft_by_amp(raft_files)
