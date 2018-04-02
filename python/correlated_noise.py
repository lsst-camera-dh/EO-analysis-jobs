import glob
from collections import namedtuple
import numpy as np
import scipy
import matplotlib.pyplot as plt
import lsst.afw.image as afw_image
import lsst.afw.math as afw_math
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest

plt.ion()
plt.rcParams['xtick.labelsize'] = 'x-small'
plt.rcParams['ytick.labelsize'] = 'x-small'


def get_overscans(infile, oscan_bbox):
    overscans = dict()
    for amp in imutils.allAmps():
        image = afw_image.ImageF(infile, imutils.dm_hdu(amp))
        overscans[amp] = image.Factory(image, oscan_bbox).getArray()
    return overscans


def get_mean_overscans(infiles, oscan_bbox):
    mean_overscans = dict()
    for amp in imutils.allAmps():
        amp_images = afw_image.vectorImageF()
        for infile in infiles:
            amp_images.push_back(afw_image.ImageF(infile, imutils.dm_hdu(amp)))
        mean_image = afw_math.statisticsStack(amp_images, afw_math.MEAN)
        mean_overscans[amp] \
            = mean_image.Factory(mean_image, oscan_bbox).getArray()
    return mean_overscans


BiasStats = namedtuple('BiasStats',
                       'noise_orig noise_corr corr_factor bias_oscan'.split())


def correlated_noise(bias_files, target=0, make_plots=False, plot_corr=True,
                     figsize=(8, 8), title=''):
    if make_plots:
        f1, ax1 = plt.subplots(4, 4, figsize=figsize)
        ax1 = {amp: subplot for amp, subplot in zip(imutils.allAmps(),
                                                    ax1.flatten())}
        f2, ax2 = plt.subplots(4, 4, figsize=figsize)
        ax2 = {amp: subplot for amp, subplot in zip(imutils.allAmps(),
                                                    ax2.flatten())}

    # Extract the target filename and omit it from the list of bias files.
    target_file = bias_files.pop(target)

    # Get the amplifier geometry and serial overscan bounding box from
    # the FITS headers.
    amp_geom = sensorTest.makeAmplifierGeometry(target_file)
    oscan_bbox = amp_geom.serial_overscan

    # Get the target frame overscans.
    bias_oscans = get_overscans(target_file, oscan_bbox)
    oscan_shape = bias_oscans[1].shape

    # Construct the mean bias overscans from the remaining files.
    mean_oscans = get_mean_overscans(bias_files, oscan_bbox)

    # Compute the mean values of the mean bias overscans.
    mean_oscan_values \
        = {amp: np.mean(oscan) for amp, oscan in mean_oscans.items()}

    # Loop over amps in target frame and compute statistics.
    bias_stats = dict()
    for amp in bias_oscans:
        print "processing amp", amp
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
                                    corr_factor, fmean1)

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
        = correlated_noise(bias_files, make_plots=True, plot_corr=False,
                           title='Run 5808D, ETU1, S00')

    bias_files = sorted(glob.glob('/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-005/6288/sflat_raft_acq/v0/37108/S00/*sflat_bias*.fits'))
    rtm005_stats, f1, f2 \
        = correlated_noise(bias_files, make_plots=True, plot_corr=False,
                           title='Run 6288, RTM-005, S00')
