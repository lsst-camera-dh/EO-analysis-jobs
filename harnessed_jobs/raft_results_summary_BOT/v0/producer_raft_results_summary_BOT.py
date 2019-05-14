#!/usr/bin/env python
"""
Producer script for BOT raft-level results summaries.
"""
import matplotlib.pyplot as plt
import siteUtils
from multiprocessor_execution import run_device_analysis_pool
from camera_components import camera_info
from focal_plane_plotting import plot_focal_plane
from bot_eo_analyses import raft_results_task, repackage_summary_files
from lsst.obs.lsst.imsim import ImsimMapper

repackage_summary_files()
run_device_analysis_pool(raft_results_task, camera_info.get_raft_names(),
                         processes=None)

#
# Focal plane heat maps
#
fp_configs = {'fe55_BOT_analysis': (('gain', 'psf_sigma'),
                                    ((0.5, 1.5), (3, 5.5)),
                                    (('e-/ADU', 'micron'))),
              'read_noise_BOT': (('read_noise',), ((0, 20),), ('e-',)),
              'bright_defects_BOT': (('bright_pixels', 'bright_columns'),
                                     (None, None), (None, None)),
              'dark_defects_BOT': (('dark_pixels', 'dark_columns'),
                                   (None, None), (None, None)),
              'dark_current_BOT': (('dark_current_95CL',), ((0, 0.3),),
                                   ('e-/pix/s',)),
              'traps_BOT': (('num_traps',), (None,), (None,)),
              'cte_BOT': (('cti_high_serial', 'cti_high_parallel',
                           'cti_low_serial', 'cti_low_parallel'),
                          ((0, 1e-5), (0, 1e-5), (0, 1e-5), (0, 1e-5)),
                          (None, None, None, None)),
              'flat_pairs_BOT': (('max_frac_dev',), ((0, 0.05),), (None,)),
              'ptc_BOT': (('ptc_gain', 'ptc_a00'), ((0.5, 1.5), (0, 5e6)),
                          ('e-/ADU', None)),
              'bright_fatter_BOT': (('bf_xcorr', 'bf_ycorr'), (None, None),
                                    (None, None))}

camera = ImsimMapper().camera

run = siteUtils.getRunNumber()
unit_id = siteUtils.getUnitId()
et_results = siteUtils.ETResults(run)
for schema_name, configs in fp_configs.items():
    for column, z_range, units in zip(*configs):
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1)
        try:
            amp_values = et_results.get_amp_data(schema_name, column)
        except KeyError:
            print("Could not find focal plane results for {}: {} in eTraveler."
                  .format(schema_name, column))
        else:
            plot_focal_plane(ax, amp_values, camera=camera, z_range=z_range)
            if units is not None:
                plt.title('Run {}, {} ({})'.format(run, column, units))
            else:
                plt.title('Run {}, {}'.format(run, column))
            outfile = '{}_{}_{}.png'.format(unit_id, run, column)
            plt.savefig(outfile)
