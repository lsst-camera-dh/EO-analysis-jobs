#!/usr/bin/env ipython
import os
import glob
import shutil
import warnings
import multiprocessing
from astropy.io import fits
import lsst.log
import lsst.eotest.sensor as sensorTest
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from bot_eo_analyses import bias_frame_task, flat_gain_stability_task, \
        make_file_prefix
from flat_gain_stability import plot_raft_by_amp


_DET_NAMES = tuple(f'R22_S{slot}' for slot in
                   '00 01 02 10 11 12 20 21 22'.split())


def find_frames(data_dir, imgtype='BIAS', testtype='FLAT',
                glob_pattern='TS_C_*/*R22_S11.fits'):
    """
    Find frames of a specified imgtype and testtype.  Return a list
    of the folder paths containing each frame.
    """
    files = sorted(glob.glob(os.path.join(data_dir, glob_pattern)))
    frame_dirs = []
    for item in files:
        with fits.open(item) as hdus:
            header = hdus[0].header
            if header['IMGTYPE'] == imgtype and header['TESTTYPE'] == testtype:
                frame_dirs.append(os.path.dirname(item))
    return frame_dirs


def make_rolloff_mask(run, det_name, template_file, outdir='.'):
    """
    Make the rolloff mask file for the given template file.
    """
    file_prefix = make_file_prefix(run, det_name)
    rolloff_mask_file = f'{file_prefix}_edge_rolloff_mask.fits'
    sensorTest.rolloff_mask(template_file, rolloff_mask_file)
    dest_dir = os.path.join(outdir, 'masks')
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, rolloff_mask_file)
    shutil.move(rolloff_mask_file, dest)
    return dest


def run_bias_frame_task(run, frame_dirs, outdir='.', processes=9,
                        det_names=_DET_NAMES, nskip=5):
    """
    Run the bias_frame_task using a multiprocessing.Pool on
    bias date for a science raft.
    """
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

    # Move data products to output directories.
    # Bias model files:
    dest_dir = os.path.join(outdir, 'bias_models')
    os.makedirs(dest_dir, exist_ok=True)
    for pattern in ('*median_bias.fits', '*pca_bias.*'):
        for src in glob.glob(pattern):
            dest = os.path.join(dest_dir, os.path.basename(src))
            shutil.move(src, dest)

    # Mask files:
    dest_dir = os.path.join(outdir, 'masks')
    os.makedirs(dest_dir, exist_ok=True)
    for pattern in ('*rolloff_mask.fits',):
        for src in glob.glob(pattern):
            dest = os.path.join(dest_dir, os.path.basename(src))
            shutil.move(src, dest)


def run_flat_gain_stability_task(run, frame_dirs, outdir='.', processes=9,
                                 det_names=_DET_NAMES):
    """
    Run the flat_gain_stability_task on flat data.
    """
    bias_dir = os.path.join(outdir, 'bias_models')
    def add_bias_dir(x):
        return os.path.join(bias_dir, x)

    mask_dir = os.path.join(outdir, 'masks')
    def add_mask_dir(x):
        return os.path.join(mask_dir, x)

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
            if os.path.isdir(bias_dir):
                bias_frame = (add_bias_dir(f'{prefix}_pca_bias.pickle'),
                              add_bias_dir(f'{prefix}_pca_bias.fits'))
                if not os.path.isfile(bias_frame[0]):
                    bias_frame = None
            else:
                bias_frame = None
            rolloff_mask_file = add_mask_dir(f'{prefix}_edge_rolloff_mask.fits')
            if os.path.isfile(rolloff_mask_file):
                mask_files = (rolloff_mask_file,)
            else:
                mask_files = (make_rolloff_mask(run, det_name, flat_files[0],
                                                outdir=outdir),)
            args = run, det_name, flat_files
            kwds = dict(mask_files=mask_files, bias_frame=bias_frame,
                        verbose=False)
            futures.append(
                pool.apply_async(flat_gain_stability_task, args, kwds))
        pool.close()
        pool.join()

        for future in futures:
            df, outfile = future.get()

    dest_dir = os.path.join(outdir, 'flat_gain_results')
    os.makedirs(dest_dir, exist_ok=True)
    for src in glob.glob(f'*_{run}_*flat_signal_sequence.pickle'):
        dest = os.path.join(dest_dir, os.path.basename(src))
        shutil.move(src, dest)


if __name__ == '__main__':
    import sys

    try:
        ts8_data_dir, outdir = sys.argv[1:3]
    except:
        print()
        print(f"usage: {os.path.basename(sys.argv[0])} <TS8 data dir> "
              "<output dir>")
        sys.exit(0)

    run = '-'.join(('TS8', os.path.basename(ts8_data_dir.rstrip('/'))))

    bias_frame_dirs = find_frames(ts8_data_dir)
    if bias_frame_dirs:
        run_bias_frame_task(run, bias_frame_dirs, outdir=outdir)

    flat_frame_dirs = find_frames(ts8_data_dir, imgtype='FLAT')
    run_flat_gain_stability_task(run, flat_frame_dirs, outdir=outdir)

    pattern = os.path.join(outdir, 'flat_gain_results',
                           f'*_{run}_flat_signal_sequence.pickle')
    raft_files = sorted(glob.glob(pattern))
    plot_raft_by_amp(raft_files)
