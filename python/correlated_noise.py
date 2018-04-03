import glob
from collections import namedtuple, defaultdict
import numpy as np
import scipy
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest

plt.ion()
plt.rcParams['xtick.labelsize'] = 'x-small'
plt.rcParams['ytick.labelsize'] = 'x-small'


def get_oscan_indices(target_file):
    "Return the pixel indices of the overscan region."
    amp_geom = sensorTest.makeAmplifierGeometry(target_file)
    bbox = amp_geom.serial_overscan
    return bbox.getMinY(), bbox.getMaxY(), bbox.getMinX(), bbox.getMaxX()


def get_overscans(infile, oscan_indices=None):
    "Return the overscan data as numpy arrays."
    if oscan_indices is None:
        y0, y1, x0, x1 = get_oscan_indices(infile)
    else:
        y0, y1, x0, x1 = oscan_indices
    ccd = fits.open(infile)
    overscans = dict()
    for amp in imutils.allAmps():
        overscans[amp] = ccd[amp].data[y0:y1, x0:x1]
    return overscans


def get_mean_overscans(infiles, oscan_indices=None):
    "Compute the mean of the overscans of the list of files"
    if oscan_indices is None:
        y0, y1, x0, x1 = get_oscan_indices(infiles[0])
    else:
        y0, y1, x0, x1 = oscan_indices
    mean_overscans = defaultdict(list)
    for infile in infiles:
        ccd = fits.open(infile)
        for amp in imutils.allAmps():
            mean_overscans[amp].append(ccd[amp].data[y0:y1, x0:x1])
    for amp, images in mean_overscans.items():
        mean_overscans[amp] = sum(images)/float(len(images))
    return mean_overscans


BiasStats = namedtuple('BiasStats',
                       'noise_orig noise_corr corr_factor bias_oscan'.split())


def correlated_noise(bias_files, target=0, make_plots=False, plot_corr=True,
                     figsize=(8, 8), title=''):
    """
    Compute the correlated noise statistics for the overscan regions
    of the list of files, optionally making plots of the distributions.
    """
    f1, f2 = None, None
    if make_plots:
        f1, ax1 = plt.subplots(4, 4, figsize=figsize)
        ax1 = {amp: subplot for amp, subplot in zip(imutils.allAmps(),
                                                    ax1.flatten())}
        f2, ax2 = plt.subplots(4, 4, figsize=figsize)
        ax2 = {amp: subplot for amp, subplot in zip(imutils.allAmps(),
                                                    ax2.flatten())}

    # Extract the target filename and omit it from the list of bias files.
    target_file = bias_files.pop(target)


    # Get the target frame overscans.
    bias_oscans = get_overscans(target_file)
    oscan_shape = bias_oscans[1].shape

    # Construct the mean bias overscans from the remaining files.
    mean_oscans = get_mean_overscans(bias_files)

    # Compute the mean values of the mean bias overscans.
    mean_oscan_values \
        = {amp: np.mean(oscan) for amp, oscan in mean_oscans.items()}

    # Loop over amps in target frame and compute statistics.
    bias_stats = dict()
    for amp in bias_oscans:
        # Loop over other amps and construct the mean image of the
        # bias-subtracted overscans.  Require included amps to have
        # (unsubtracted) overscans with 4 < stdev < 25 rms ADU.
        reduced_mean_oscan = np.zeros(oscan_shape)
        num_oscan = 0
        for oamp, oscan in bias_oscans.items():
            if oamp == amp or not (4. < np.std(oscan) < 25):
                continue
            reduced_mean_oscan += (oscan - mean_oscans[oamp])
            num_oscan += 1
        reduced_mean_oscan -= np.mean(reduced_mean_oscan)
        reduced_mean_oscan /= num_oscan

        fdata1 = bias_oscans[amp] - mean_oscans[amp]
        fmean1 = np.mean(fdata1)
        fdata1 -= fmean1
        dmat = np.vstack((reduced_mean_oscan.flatten(), fdata1.flatten()))
        covmat = scipy.cov(dmat, rowvar=True)
        corr_factor = covmat[0, 1]/covmat[0, 0]
        fdiff = fdata1 - corr_factor*reduced_mean_oscan
        bias_stats[amp] = BiasStats(np.sqrt(covmat[1, 1]), np.std(fdiff),
                                    corr_factor,
                                    np.mean(bias_oscans[amp]))
                                    #fmean1)

        if make_plots:
            f1.suptitle(title)
            f2.suptitle(title)
            ax1[amp].hist2d(reduced_mean_oscan.flatten(), fdata1.flatten(),
                            bins=(100, 100), range=((-50, 50), (-50, 50)))
            if plot_corr:
                ax2[amp].hist(fdiff.flatten(), bins=100, range=(-50, 50),
                              histtype='step')
            else:
                ax2[amp].hist(fdata1.flatten(), bins=100, range=(-50, 50),
                              histtype='step')

    return bias_stats, f1, f2


if __name__ == '__main__':
    bias_files = sorted(glob.glob('/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive-test/LCA-11021_RTM/LCA-11021_RTM-002_ETU1/5808D/sflat_raft_acq/v0/38318/S00/ITL-3800C-023-Dev_sflat_bias_*.fits'))
    etu1_stats, f1, f2 \
        = correlated_noise(bias_files, make_plots=False, plot_corr=False,
                           title='Run 5808D, ETU1, S00')
    for amp, stats in etu1_stats.items():
        print amp, stats

    bias_files = sorted(glob.glob('/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-005/6288/sflat_raft_acq/v0/37108/S00/*sflat_bias*.fits'))
    rtm005_stats, f1, f2 \
        = correlated_noise(bias_files, make_plots=True, plot_corr=False,
                           title='Run 6288, RTM-005, S00')
    for amp, stats in rtm005_stats.items():
        print amp, stats
