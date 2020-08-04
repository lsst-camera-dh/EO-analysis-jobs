import os
from collections import defaultdict, namedtuple
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from eTraveler.clientAPI.connection import Connection
import lsst.eotest.image_utils as imutils
from camera_components import camera_info
import siteUtils
from focal_plane_plotting import plot_focal_plane
plt.ion()

# Snap-shot of good TS8 runs listed at
# https://confluence.slac.stanford.edu/display/LSSTCAM/List+of+Good+Runs
ts8_good_runs = {'LCA-11021_RTM-018': '12120',
                 'LCA-11021_RTM-010': '12139',
                 'LCA-11021_RTM-021': '12086',
                 'LCA-11021_RTM-016': '12027',
                 'LCA-11021_RTM-015': '12002',
                 'LCA-11021_RTM-004': '11977',
                 'LCA-11021_RTM-008': '11952',
                 'LCA-11021_RTM-007': '11903',
                 'LCA-11021_RTM-005': '11852',
                 'LCA-11021_RTM-019': '11808',
                 'LCA-11021_RTM-006': '11746',
                 'LCA-11021_RTM-022': '11671',
                 'LCA-11021_RTM-009': '11415',
                 'LCA-11021_RTM-024': '11351',
                 #'LCA-11021_RTM-017': '11166',
                 'LCA-11021_RTM-017': '11188',
                 'LCA-11021_RTM-012': '11063',
                 'LCA-11021_RTM-013': '10982',
                 'LCA-11021_RTM-014': '10928',
                 'LCA-11021_RTM-023': '10517',
                 'LCA-11021_RTM-020': '10669',
                 'LCA-11021_RTM-025': '10722',
                 'LCA-11021_RTM-011': '10861'}

ts8_gain_scale_factors = {'R41': 1.5, 'R42': 1.5, 'R10': 1.5}

RaftInfo = namedtuple('RaftInfo', ['lsst_id', 'good_run'])


def get_raft_slot_info():
    """
    Return a dict, keyed by raft slot (e.g., 'R01') of tuples containing
    raft LSST ID (e.g., 'LCA-11021_RTM-017') and the "good" TS8 run number
    from
    https://confluence.slac.stanford.edu/display/LSSTCAM/List+of+Good+Runs
    """
    cryostat_id = 'LCA-10134_Cryostat-0001'
    cryo_htype = 'LCA-10134_Cryostat'
    conn = Connection('jchiang', 'Prod', prodServer=True)
    resp = conn.getHardwareHierarchy(experimentSN=cryostat_id, htype=cryo_htype)

    raft_htype = 'LCA-11021_RTM'
    raft_slot_info = dict()
    for item in resp:
        if item['child_hardwareTypeName'] == raft_htype:
            raft_slot = f'R{item["slotName"][-2:]}'
            lsst_id = item['child_experimentSN']
            raft_slot_info[raft_slot] \
                = RaftInfo(lsst_id, ts8_good_runs[lsst_id])
    return raft_slot_info


def get_ts8_amp_data(schema_name, field_name):
    """
    Get the amp values for the specified schema_name and field_name
    (e.g., 'fe55_raft_analysis' and 'gain') for the good TS8 runs.
    """
    raft_info = get_raft_slot_info()
    amp_data = defaultdict(dict)
    for raft_slot, (_, run) in raft_info.items():
        df = siteUtils.ETResults(run=run)[schema_name]
        for _, row in df.iterrows():
            det_name = f'{raft_slot}_{row.slot}'
            channel = 'C' + imutils.channelIds[row.amp]
            amp_data[det_name][channel] = row[field_name]
            if raft_slot in ts8_gain_scale_factors:
                amp_data[det_name][channel] *= ts8_gain_scale_factors[raft_slot]
    return amp_data


def get_bot_amp_data(run, schema_name, field_name):
    """
    Get the amp values for the specified run, schema_name, and field_name
    for BOT data.
    """
    results = siteUtils.ETResults(run=run)
    return results.get_amp_data(schema_name, field_name)


def adjust_gains(amp_gains, gain_max=1.5):
    """Set some corner raft gains by hand."""
    my_amp_gains = dict()
    my_amp_gains['R04_SG0'] = dict(C13=1.25, C14=1.25, C16=1.25, C17=1.25)
    my_amp_gains['R40_SG1'] = dict(C10=1.18, C11=1.18, C17=1.18)
    my_amp_gains['R44_SG0'] = dict(C05=1.17)
    for detname, gains in my_amp_gains.items():
        for channel, gain in gains.items():
            if amp_gains[detname][channel] > gain_max:
                amp_gains[detname][channel] = gain
    return amp_gains


if __name__ == '__main__':
    amp_gain_file = 'bot_ts8_amp_gains.json'
    if not os.path.isfile(amp_gain_file):
        amp_data = dict(
            ts8_fe55_gains=get_ts8_amp_data('fe55_raft_analysis', 'gain'),
            ts8_ptc_gains=get_ts8_amp_data('ptc_raft', 'ptc_gain'),
            bot_fe55_gains=get_bot_amp_data('6801D', 'fe55_BOT_analysis',
                                            'gain'))

        with open(amp_gain_file, 'w') as output:
            json.dump(amp_data, output)
    else:
        with open(amp_gain_file, 'r') as fd:
            amp_data = json.load(fd)

    gain_min = 0.5
    gain_max = 1.6
    amp_gains = defaultdict(dict)
    for det_name, gain_values in amp_data['ts8_fe55_gains'].items():
        for channel, fe55_gain in gain_values.items():
            ptc_gain = amp_data['ts8_ptc_gains'][det_name][channel]
            if gain_min < fe55_gain < gain_max:
                amp_gains[det_name][channel] = fe55_gain
            else:
                amp_gains[det_name][channel] = ptc_gain
        if det_name in amp_data['bot_fe55_gains']:
            bot_gains = amp_data['bot_fe55_gains'][det_name]
            for channel, bot_gain in bot_gains.items():
                if gain_min < bot_gain < gain_max:
                    amp_gains[det_name][channel] = bot_gain
    for det_name, corner_raft_gains in amp_data['bot_fe55_gains'].items():
        if det_name not in amp_gains:
            amp_gains[det_name].update(corner_raft_gains)

    # Adjust some corner raft gains by hand:
    amp_gains = adjust_gains(amp_gains)

    curated_amp_gains = 'curated_amp_gains.json'
    with open(curated_amp_gains, 'w') as output:
        json.dump(amp_gains, output)

    for gain_source, my_gains in amp_data.items():
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111)
        plot_focal_plane(ax, my_gains, camera=camera_info.camera_object,
                         z_range=(gain_min, gain_max))
        plt.title(gain_source)
        plt.savefig(f'{gain_source}.png')

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    plot_focal_plane(ax, amp_gains, camera=camera_info.camera_object,
                     z_range=(gain_min, gain_max))
    plt.title('combined ts8 fe55, ts8 ptc, bot_fe55 gains')
    plt.savefig(f'amp_gains.png')

