"""
Module to compute corrections to the photodiode integrals based on
mismatches between signal at boundaries between filter combinations in
a flat pair sequence.
"""
import os
import glob
from collections import defaultdict
import numpy as np
from astropy.io import fits
import pandas as pd


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


def flat_metadata(flat_dir):
    """
    Extract frame metadata from the symlink name and physical path of
    the folder containing the flat files.
    """
    # Get the target signal and filter combination from the symlink name.
    tokens = os.path.basename(flat_dir)[len('flat_'):-len('_flat0_000')]\
                    .split('_')
    signal = float(tokens[-1])
    filt = '_'.join(tokens[:-1])

    # Get the dayobs and seqnum from the physical path.
    tokens = os.path.realpath(flat_dir).split('_')
    seqnum = int(tokens[-1])
    dayobs = int(tokens[-2])

    return signal, filt, seqnum, dayobs


def flat_sequence_metadata(data_dir):
    """
    Return a data frame with the metadata from the flat folders in
    the specified data directory.
    """
    data = defaultdict(list)
    flats = glob.glob(os.path.join(data_dir, 'flat*flat?_*'))
    for flat in flats:
        signal, filt, seqnum, dayobs = flat_metadata(flat)
        data['signal'].append(signal)
        data['filt'].append(filt)
        data['seqnum'].append(seqnum)
        data['dayobs'].append(dayobs)
    return pd.DataFrame(data=data)


def apply_filter_corrections(detresp_file, data_dir, amp=1, max_signal=3e4):
    """
    Compute and apply per filter corrections to photodiode integrals (flux)
    in the detresp_file.
    """
    # Use the specified amp for the signal level to use for calibration.
    xcol, ycol = 'flux', f'AMP{amp:02d}_SIGNAL'

    # Extract the metadata from the flat pair sequence in order
    # to associate the filter combination for each frame.
    md = flat_sequence_metadata(data_dir)

    # Index each filter combo by DAYOBS and SEQNUM.
    filters = {(d, s): f for d, s, f in
               zip(md['dayobs'], md['seqnum'], md['filt'])}

    # Read the detector response file into a data frame.
    with fits.open(detresp_file) as hdus:
        data = {col.name: [_ for _ in hdus[1].data[col.name]]
                for col in hdus[1].data.columns}
        df0 = pd.DataFrame(data=data)

    # Add the filter info.
    df0['filt'] = [filters.get((dayobs, seqnum), 'None')
                   for dayobs, seqnum in zip(df0['DAYOBS'], df0['SEQNUM'])]
    df = pd.DataFrame(df0.query(f'filt != "None" and {ycol} < {max_signal}'))

    # Compute the correction factors.
    corr_factors, pd_corrections \
        = pd_filter_corrections(df[xcol].to_numpy(), df[ycol].to_numpy(),
                                df['filt'].to_numpy())

    # Apply the corrections.
    df[xcol] *= corr_factors

    # Return a tuple of the corrected detector response data and the
    # corrections for each filter combination.
    return df, pd_corrections
