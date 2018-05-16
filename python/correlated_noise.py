import glob
from collections import namedtuple, defaultdict
import itertools
import numpy as np
import scipy
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import camera_components

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

    Parameters
    ----------
    bias_files: list
        List of bias files to analyze.  This list must have at least as many
        files as the target file index + 1.
    target: int
        Bias frame to compare to the mean biases constructed from the
        remaining files.
    make_plots: bool [False]
        Flag to determine if the png plots will be generated.
    plot_corr: bool [True]
        Flag to plot the histograms of correlation-corrected pixel
        values.  If False, then plot histograms of the uncorrected pixel
        values.
    figsize: tuple [(8, 8)]
        Figure size (in inches) of 4x4 grid of correlation plots.
    title: str ['']
        Title of 4x4 grid.

    Returns
    -------
    (dict, figure, figure):  tuple of results and matplotlib
        figures.  The first item is a dict of BiasStats objects,

        BiasStats = namedtuple('BiasStats', \
                    'noise_orig noise_corr corr_factor bias_oscan'.split())

        that contain the results for each amplifier.
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
            label = 'amp %i, cov/var = %.2f' \
                    % (amp, bias_stats[amp].corr_factor)
            ax1[amp].text(-40, 40, label, fontsize=6, color='w',
                          fontweight='bold')
            if plot_corr:
                ax2[amp].hist(fdiff.flatten(), bins=100, range=(-50, 50),
                              histtype='step')
            else:
                ax2[amp].hist(fdata1.flatten(), bins=100, range=(-50, 50),
                              histtype='step')

    return bias_stats, f1, f2

def raft_level_oscan_correlations(bias_files, buffer=10, title='',
                                  vmin=0, vmax=0.5):
    """
    Compute the correlation coefficients between the overscan pixels
    of the 144 amplifiers in raft.

    Parameters
    ----------
    bias_files: dict
        Dictionary of bias image files, indexed by sensor slot id.
    buffer: int [10]
        Buffer region around perimeter of serial overscan region to
        avoid when computing the correlation coefficients.
    title: str ['']
        Plot title.
    vmin: float [0]
        Minimum pixel value for color scale.
    vmax: float [0.5]
        Maximum pixel value for color scale.

    Returns
    -------
    (matplotlib.figure.Figure, np.array): The figure containing the plot and
        the numpy array containing the correlation coefficients.
    """
    slots = 'S00 S01 S02 S10 S11 S12 S20 S21 S22'.split()
    bbox = None
    overscans = []
    for slot in slots:
        ccd = sensorTest.MaskedCCD(bias_files[slot])
        if bbox is None:
            bbox = ccd.amp_geom.serial_overscan
            bbox.grow(-buffer)
        for amp in ccd:
            image = ccd[amp].getImage()
            overscans.append(image.Factory(image, bbox).getArray())
    namps = len(overscans)
    data = np.array([np.corrcoef(overscans[i[0]].ravel(),
                                 overscans[i[1]].ravel())[0, 1]
                     for i in itertools.product(range(namps), range(namps))])
    data = data.reshape((namps, namps))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize='medium')
    image = ax.imshow(data, interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.colorbar(image)
    set_ticks(ax, slots, amps=16)
    return fig, data

def set_ticks(ax, slots, amps=16):
    """Set the tick labels, centering the slot names between amps 1 and 16."""
    major_locs = [i*amps - 0.5 for i in range(len(slots) + 1)]
    minor_locs = [amps//2 + i*amps for i in range(len(slots))]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_tick_params(which='minor', length=0)
        axis.set_major_locator(ticker.FixedLocator(major_locs))
        axis.set_major_formatter(ticker.FixedFormatter(['']*len(major_locs)))
        axis.set_minor_locator(ticker.FixedLocator(minor_locs))
        axis.set_minor_formatter(ticker.FixedFormatter(slots))

if __name__ == '__main__':
    plt.ion()
    bias_files = sorted(glob.glob('/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive-test/LCA-11021_RTM/LCA-11021_RTM-002_ETU1/5808D/sflat_raft_acq/v0/38318/S00/ITL-3800C-023-Dev_sflat_bias_*.fits'))
    etu1_stats, f1, f2 \
        = correlated_noise(bias_files, make_plots=True, plot_corr=False,
                           title='Run 5808D, ETU1, S00')
    for amp, stats in etu1_stats.items():
        print amp, stats
    plt.figure(f1.number)
    plt.savefig('ETU1_S00_noise_corr.png')

    bias_files = sorted(glob.glob('/gpfs/slac/lsst/fs1/g/data/jobHarness/jh_archive/LCA-11021_RTM/LCA-11021_RTM-005/6288/sflat_raft_acq/v0/37108/S00/*sflat_bias*.fits'))
    rtm005_stats, f1, f2 \
        = correlated_noise(bias_files, make_plots=True, plot_corr=False,
                           title='Run 6288, RTM-005, S00')
    for amp, stats in rtm005_stats.items():
        print amp, stats
    plt.figure(f1.number)
    plt.savefig('RTM-005_S00_noise_corr.png')
