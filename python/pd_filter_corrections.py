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

    # Add the filter info
    df0['filt'] = [filters.get((dayobs, seqnum), 'None')
                   for dayobs, seqnum in zip(df0['DAYOBS'], df0['SEQNUM'])]
    df = pd.DataFrame(df0.query(f'filt != "None" and {ycol} < {max_signal}'))

    # Sort filters by flux, in ascending order.
    filters = sorted(list(set(df['filt'])),
                     key=lambda x: np.median(df.query(f'filt=="{x}"')[xcol]))

    # Fit the linear scale factor independently for the data associated
    # with each filter combo.
    y0 = dict()
    for filt in filters:
        my_df = df.query(f'filt=="{filt}"')
        y0[filt] = sum(my_df[ycol])/sum(my_df[xcol])

    # Compute corrections relative to the filter combination of the
    # highest signal data.
    pd_corrections = {filt: y0[filt]/y0[filters[-1]] for filt in filters}

    # Apply the corrections.
    df[xcol] *= np.array([pd_corrections.get(filt, 1) for filt in df['filt']])

    # Return a tuple of the corrected detector response data and the
    # corrections for each filter combination.
    return df, pd_corrections
