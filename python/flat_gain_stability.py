"""
Module to analyze flat gain stability data.
"""
import os
import glob
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import siteUtils


__all__ = ['plot_raft', 'plot_all_rafts']


def plot_raft(raft_files, cut=None, divide_by_flux=True, y_range=None,
              colors=None):
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    raft = os.path.basename(raft_files[0])[:len('R22')]
    mjd_min = None
    for item, color in zip(raft_files, colors):
        sensor = os.path.basename(item)[len('R22_'):len('R22_S11')]
        raw_df = pd.read_pickle(item)
        df = raw_df.query(cut) if cut is not None else raw_df
        if mjd_min is None:
            mjd_min = np.int(min(df['mjd']))
            flux = df['flux']/np.mean(df['flux'])
            time = (24*(df['mjd'] - mjd_min)).to_list()
        x, y = [], []
        for amp in range(1, 17):
            try:
                signal = df[f'amp{amp:02d}']
            except KeyError:
                continue
            if divide_by_flux:
                signal /= flux
            signal /= np.mean(signal)
            x.extend(time)
            y.extend(signal.to_list())
        plt.scatter(x, y, s=1, color=color, label=sensor)
    plt.xlabel(f'(MJD - {mjd_min})*24 (hours)')
    plt.ylabel('counts/mean(counts)')
    plt.legend(fontsize='x-small')
    plt.title(raft, fontsize='small')
    if y_range is not None:
        plt.ylim(*y_range)


def plot_all_rafts(run, results_dir='.', cut=None, y_range=(0.998, 1.002),
                   divide_by_flux=True, figsize=(18, 18)):
    """
    Plot the flat gain stability curve for all 25 rafts in a 5x5 grid.
    """
    raft_files = defaultdict(list)
    files = sorted(glob.glob(os.path.join(results_dir,
                                          '*flat_signal_sequence.pickle')))
    if not files:
        return
    for item in files:
        raft_name = os.path.basename(item)[:len('R22')]
        raft_files[raft_name].append(item)

    figure = plt.figure(figsize=figsize)
    rafts = sorted(list(raft_files.keys()))
    for i, raft in enumerate(rafts):
        ax = figure.add_subplot(5, 5, i+1)
        plot_raft(raft_files[raft], cut=cut, divide_by_flux=divide_by_flux,
                  y_range=y_range)
    plt.tight_layout()
    plt.suptitle(f'Run {run}')
    unit_id = siteUtils.getUnitId()
    plt.savefig(f'{unit_id}_{run}_flat_gain_stability.png')
