"""
Module to make Fe55 gain stability plots.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import siteUtils
import lsst.eotest.image_utils as imutils
from bot_eo_analyses import raft_ccd_order, append_acq_run

__all__ = ['plot_raft_fe55_gains_by_amp', 'plot_raft_fe55_gains_by_ccd',
           'plot_all_raft_fe55_gains']


channel = {i: f'C{_}' for i, _ in imutils.channelIds.items()}


def plot_raft_fe55_gains_by_amp(raft_files, figsize=(12, 12), y_range=None,
                                outfile=None):
    """
    Plot Fe55 gains for all amps for each CCD in a raft.

    Parameters
    ----------
    raft_files: list
        A list of pickle files containing pandas data frames with
        the gains for each amp in a CCD as a function of MJD.
        The filenames are assumed to be of the form
        '<raft>_<sensor>_<run>_gain_sequence.pickle`, so that the
        detector name (e.g., 'R22_S11'), raft, and run number can
        be extracted.
    figsize: (float, float) [(12, 12)]
        Size of the figure for each raft.
    y_range: (float, float) [None]
        Plotting limts of the y-axis for each plot.  If None, then the
        y-limits will be scaled by matplotlib to the data.
    outfile: str [None]
        Filename of output png file. If None, then a default name of
        f'{raft}_{run}_fe55_gain_stability.png' will be used.
    """
    sorted_raft_files = raft_ccd_order.sorted(raft_files)
    plt.figure(figsize=figsize)
    for i, item in enumerate(sorted_raft_files, 1):
        plt.subplot(3, 3, i)
        df = pd.read_pickle(item)
        mjd0 = int(min(df['mjd']))
        det_name = os.path.basename(item)[:len('R22_S11')]
        amps = sorted(list(set(df['amp'])))
        for amp in amps:
            my_df = df.query(f'amp == {amp}')
            gains = my_df['gain'].to_numpy()
            frac = gains/np.mean(gains)
            plt.scatter(24*(my_df['mjd'] - mjd0), frac, s=2,
                        label=f'{channel[amp]}')
        plt.legend(fontsize='x-small', ncol=2)
        if y_range is not None:
            plt.ylim(*y_range)
        plt.title(det_name, fontsize='x-small')
        plt.xlabel(f'(MJD - {mjd0:d})*24 (hours)')
        plt.ylabel('gain/mean(gain)')
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    tokens = os.path.basename(raft_files[0]).split('_')
    run = tokens[2]
    raft = tokens[0]
    suptitle = append_acq_run(f'Fe55 gain stability, {raft}, Run {run}')
    plt.suptitle(suptitle)
    if outfile is None:
        outfile = f'{raft}_{run}_fe55_gain_stability.png'
    plt.savefig(outfile)


def plot_raft_fe55_gains_by_ccd(raft_files, colors=None, y_range=None):
    """
    Plot Fe55 gains for each ccd in a raft, aggregating over amps.
    """
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    raft = os.path.basename(raft_files[0])[:len('R22')]
    mjd0 = None
    sorted_raft_files = raft_ccd_order.sorted(raft_files)
    for item, color in zip(sorted_raft_files, colors):
        sensor = os.path.basename(item)[len('R22_'):len('R22_S11')]
        df = pd.read_pickle(item)
        if mjd0 is None:
            mjd0 = int(min(df['mjd']))
        x, y = [], []
        for amp in range(1, 17):
            my_df = df.query(f'amp == {amp}')
            time = 24*(my_df['mjd'] - mjd0)
            gain = my_df['gain']/np.mean(my_df['gain'])
            x.extend(time)
            y.extend(gain)
        plt.scatter(x, y, s=1, color=color, label=sensor)
    plt.xlabel(f'(MJD - {mjd0})*24 (hours)')
    plt.ylabel('gain/mean(gain)')
    plt.legend(fontsize='x-small', ncol=2)
    plt.title(raft, fontsize='small')
    if y_range is not None:
        plt.ylim(*y_range)


def plot_all_raft_fe55_gains(raft_files, figsize=(18, 18), y_range=None):
    """
    Plot the flat gain stability curve for all 25 rafts in a 5x5 grid.
    """
    figure = plt.figure(figsize=figsize)
    rafts = sorted(list(raft_files.keys()))
    for i, raft in enumerate(rafts, 1):
        figure.add_subplot(5, 5, i)
        plot_raft_fe55_gains_by_ccd(raft_files[raft], y_range=y_range)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    run = os.path.basename(raft_files[rafts[0]][0]).split('_')[2]
    suptitle = append_acq_run(f'Fe55 gain stability, Run {run}')
    plt.suptitle(suptitle)
    unit_id = siteUtils.getUnitId()
    plt.savefig(f'{unit_id}_{run}_fe55_gain_stability.png')
