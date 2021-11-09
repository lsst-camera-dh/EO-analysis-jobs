"""
Module to compute corrections to the photodiode integrals based on
mismatches between signal at boundaries between filter combinations in
a flat pair sequence.
"""
import os
from collections import defaultdict
import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from astropy.io import fits


def pd_filter_corrections(flux, Ne, filters, Ne_max=3e4):
    """
    Function to compute filter-specific corrections to the photodiode
    integrals, bootstrapping from the data itself.

    Parameters
    ----------
    flux: np.array
        Photodiode integral values for each exposure in the sequence.
    Ne: np.array
        Measured signal level from the CCDs.  This will nominally be the
        e-/pixel level for an amp in the CCD being analyzed.
    filters: np.array
        Array of filter combinations, e.g., 'empty_SDSSi', 'ND_OD0.1_SDSSi'.
    Ne_max: float [3e4]
        Maximum signal level (e-/pixel) to consider for fitting linear
        scale factors for each filter combo. This should avoid low enough
        to avoid rollover from B/F or nonlinearity at high signal levels.

    Returns
    -------
    (np.array, dict)  The first element is a np.array of the correction
    factors to apply directly to the flux array, and second element is a
    dictionary of correction factors, keyed by filter combination.
    """
    # Determine filter order by sorting on signal level.
    filter_list = []
    for filt in filters[np.argsort(Ne)]:
        if filt not in filter_list:
            filter_list.append(filt)

    # Fit the linear scale factor, y0, assuming the model Ne = y0*flux,
    # independently for each filter under the requirement Ne < Ne_max.
    y0 = dict()
    for filt in filter_list:
        index = np.where((filters == filt) & (Ne < Ne_max))
        y0[filt] = sum(Ne[index])/sum(flux[index])

    # Compute the corrections relative to the filter combination of the
    # highest signal data.
    pd_corrections = {filt: y0[filt]/y0[filter_list[-1]]
                      for filt in filter_list}

    # Return a numpy array of the correction factors and a dictionary
    # of those corrections, keyed by filter combination.
    return (np.array([pd_corrections[filt] for filt in filters]),
            pd_corrections)


def apply_corrections(fluxes, filters, pd_corrections):
    """Apply pd correctionst to the fluxes by fitler."""
    return np.array([flux*pd_corrections.get(filt, 1)
                     for flux, filt in zip(fluxes, filters)])


def plot_pd_corrections(det_resp_files, x_range=(0.975, 1.015),
                        png_file=None, pd_corr_file=None,
                        pd_corrections=None):
    """
    Plot distributions of photodiode filter corrections derived from
    a list of detector response files, which are all assumed to
    be from the same analysis run.
    """
    my_pd_corrections = defaultdict(list)
    for item in det_resp_files:
        #print(item)
        with fits.open(item) as det_resp:
            filters = det_resp[1].data['filter']
            index = np.where(filters != '')
            filters = filters[index]
            flux = det_resp[1].data['flux'][index]
            if pd_corrections is not None:
                flux = apply_corrections(flux, filters, pd_corrections)
            if '_SW' in os.path.basename(item):
                amps = range(1, 9)
            else:
                amps = range(1, 17)
            for amp in amps:
                Ne = det_resp[1].data[f'AMP{amp:02d}_SIGNAL'][index]
                try:
                    _, pd_corrs = pd_filter_corrections(flux, Ne, filters)
                except ZeroDivisionError:
                    print(os.path.basename(item), amp)
                    continue
                for filt, value in pd_corrs.items():
                    my_pd_corrections[filt].append(value)
    run = os.path.basename(det_resp_files[0]).split('_')[2]
    if pd_corrections is None:
        if pd_corr_file is None:
            pd_corr_file = f'pd_corrections_{run}.pickle'
        with open(pd_corr_file, 'wb') as fd:
            pickle.dump(my_pd_corrections, fd)

    plt.figure()
    bins = 40
    for filt, values in my_pd_corrections.items():
        plt.hist(values, alpha=0.5, bins=bins, label=filt, range=x_range)
        try:
            est_bw = (x_range[1] - x_range[0])/bins*50
            kernel = scipy.stats.gaussian_kde(values, bw_method=est_bw)
            xvals = np.linspace(x_range[0], x_range[1], 1000)
            yvals = kernel(xvals)
            x_mode = xvals[np.where(yvals == max(yvals))][0]
            plt.plot(xvals, yvals, linestyle='--', color='black', alpha=0.5)
        except np.linalg.LinAlgError:
            x_mode = np.median(values)
        print(filt, np.mean(values), np.median(values), x_mode)
        plt.axvline(x_mode, linestyle=':', color='black', alpha=0.5)
    plt.xlabel('pd correction factor')
    plt.ylabel('entries / bin')
    plt.title(f'Run {run}')
    plt.legend(fontsize='x-small')
    if png_file is None:
        png_file = f'pd_corrections_{run}.png'
    plt.savefig(png_file)
