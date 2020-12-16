import copy
import itertools
import numpy as np
import matplotlib.pyplot as plt
import astropy.visualization as viz
from astropy.visualization.mpl_normalize import ImageNormalize
import lsst.eotest.sensor as sensorTest
from .correlated_noise import set_ticks

def diff_image_arrays(flat1_file, flat2_file, bias_frame, buffer=10):
    """
    Compute the difference images for each amp, using the Astier
    weighting scheme to account for somewhat different exposure times,
    and return a dict of the image arrays, keyed by amp.
    """
    image_arrays = dict()
    ccd1 = sensorTest.MaskedCCD(flat1_file, bias_frame=bias_frame)
    ccd2 = sensorTest.MaskedCCD(flat2_file, bias_frame=bias_frame)
    imaging_bbox = ccd1.amp_geom.imaging
    imaging_bbox.grow(-buffer)
    for amp in ccd1:
        image1 = ccd1.unbiased_and_trimmed_image(amp, imaging=imaging_bbox)
        image2 = ccd2.unbiased_and_trimmed_image(amp, imaging=imaging_bbox)
        mean1 = afwMath.makeStatistics(image1, afwMath.MEAN, ccd1.stat_ctrl)\
                       .getValue()
        mean2 = afwMath.makeStatistics(image2, afwMath.MEAN, ccd1.stat_ctrl)\
                       .getValue()
        fmean = (mean1 + mean2)/2.
        image1 *= mean2/fmean
        image2 *= mean1/fmean
        image_arrays[amp] = copy.deepcopy((image1 - image2).getImage()
                                          .getArray())
    return image_arrays


def raft_level_signal_correlations(flat1_files, flat2_files, bias_frames,
                                   buffer=10, title='', vrange=None,
                                   stretch=viz.LinearStretch, figsize=(8, 8)):

    """
    Compute the correlation coefficients between the imaging section pixels
    for the difference images from a flat pair for the 144 amplifiers in raft.

    Parameters
    ----------
    flat1_files: dict
        Dictionary of flat1 image files, indexed by sensor slot id.
        These should be from the same flat pair frame as the flat2_files.
    flat2_files: dict
        Dictionary of flat2 image files, indexed by sensor slot id.
    bias_frames: dict
        Dictionary of super bias frames, indexed by sensor slot id.
    buffer: int [10]
        Buffer region around perimeter of serial overscan region to
        avoid when computing the correlation coefficients.
    title: str ['']
        Plot title.
    vrange: (float, float) [None]
        Minimum and maximum values for color scale range. If None, then
        the range of the central 98th percentile of the absolute value
        of the data is used.
    stretch: astropy.visualization.BaseStretch [LinearStretch]
        Stretch to use for the color scale.

    Returns
    -------
    (matplotlib.figure.Figure, np.array): The figure containing the plot and
        the numpy array containing the correlation coefficients.
    """
    slots = 'S00 S01 S02 S10 S11 S12 S20 S21 S22'.split()
    segments = []

    ccd0 = sensorTest.MaskedCCD(list(flat1_files.values())[0])
    bbox = ccd0.amp_geom.imaging
    bbox.grow(-buffer)

    for slot in slots:
        if slot not in flat1_files:
            for amp in ccd0:
                segments.append(np.zeros((bbox.getHeight(), bbox.getWidth())))
        else:
            imarrs = diff_image_arrays(flat1_files[slot], flat2_files[slot],
                                       bias_frame=bias_frames[slot],
                                       buffer=buffer)
            for amp in imarrs:
                segments.append(imarrs[amp])
    namps = len(segments)
    data = np.array([np.corrcoef(segments[i[0]].ravel(),
                                 segments[i[1]].ravel())[0, 1]
                     for i in itertools.product(range(namps), range(namps))])
    data = data.reshape((namps, namps))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize='medium')

    interval = viz.PercentileInterval(98.)
    if vrange is None:
        vrange = interval.get_limits(np.abs(data.ravel()))
    norm = ImageNormalize(vmin=vrange[0], vmax=vrange[1], stretch=stretch())
    image = ax.imshow(data, interpolation='none', norm=norm)
    plt.colorbar(image)

    set_ticks(ax, slots, amps=16)

    return fig, data
