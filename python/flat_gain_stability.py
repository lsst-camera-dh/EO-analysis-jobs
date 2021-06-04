"""
Module to analyze flat gain stability data.
"""
import os
import glob
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import lsst.afw.math as afwMath
import lsst.eotest.image_utils as imutils
import siteUtils
from bot_eo_analyses import raft_ccd_order, append_acq_run

__all__ = ['plot_raft_by_amp', 'plot_raft', 'plot_all_rafts']


channel = {i: f'C{_}' for i, _ in imutils.channelIds.items()}


def discard_outliers(signal, nsig=5):
    stats = afwMath.makeStatistics(signal, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    mean, stdev = (stats.getValue(afwMath.MEANCLIP),
                   stats.getValue(afwMath.STDEVCLIP))
    return np.where(np.abs(signal - mean) < nsig*stdev)


def plot_raft_by_amp(raft_files, cut=None, divide_by_flux=True,
                     figsize=(12, 12), y_range=None, outfile=None):
    """
    Plot flat gain stability for each amp individually.
    """
    sorted_raft_files = raft_ccd_order.sorted(raft_files)
    plt.figure(figsize=figsize)
    for i, item in enumerate(sorted_raft_files, 1):
        plt.subplot(3, 3, i)
        df = pd.read_pickle(item)
        if cut:
            df = df.query(cut)
        flux = (df['flux']/np.mean(df['flux'])).to_numpy()
        mjd0 = int(min(df['mjd']))
        det_name = os.path.basename(item)[:len('R22_S11')]
        for amp in range(1, 17):
            try:
                signal = df[f'amp{amp:02d}'].to_numpy()
            except KeyError:
                continue
            if divide_by_flux:
                signal /= flux
            index = discard_outliers(signal)
            signal /= np.mean(signal[index])
            times = 24*(df['mjd'].to_numpy() - mjd0)[index]
            plt.scatter(times, signal[index], s=2, label=channel[amp])
        plt.legend(fontsize='x-small', ncol=2)
        if y_range is not None:
            plt.ylim(*y_range)
        plt.title(det_name, fontsize='x-small')
        plt.xlabel(f'(MJD - {mjd0})*24 (hours)')
        plt.ylabel('counts/mean(counts)')
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    tokens = os.path.basename(raft_files[0]).split('_')
    run = tokens[2]
    raft = tokens[0]
    suptitle = append_acq_run(f'Flat gain stability, {raft}, Run {run}')
    plt.suptitle(suptitle)
    if outfile is None:
        outfile = f'{raft}_{run}_flat_gain_stability.png'
    plt.savefig(outfile)


def plot_raft(raft_files, cut=None, divide_by_flux=True, y_range=None,
              colors=None):
    """
    Plot flat gain stability for all ccds in a raft, aggregating over amps.
    """
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    raft = os.path.basename(raft_files[0])[:len('R22')]
    mjd_min = None
    sorted_raft_files = raft_ccd_order.sorted(raft_files)
    for item, color in zip(sorted_raft_files, colors):
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
    plt.legend(fontsize='x-small', ncol=2)
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
        figure.add_subplot(5, 5, i+1)
        plot_raft(raft_files[raft], cut=cut, divide_by_flux=divide_by_flux,
                  y_range=y_range)
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    suptitle = append_acq_run(f'Flat gain stability, Run {run}')
    plt.suptitle(suptitle)
    unit_id = siteUtils.getUnitId()
    plt.savefig(f'{unit_id}_{run}_flat_gain_stability.png')
