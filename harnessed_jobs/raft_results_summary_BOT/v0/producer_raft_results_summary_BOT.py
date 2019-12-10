#!/usr/bin/env ipython
"""
Producer script for BOT raft-level results summaries.
"""
import os
import matplotlib.pyplot as plt
import siteUtils
from camera_components import camera_info
from focal_plane_plotting import plot_focal_plane, hist_amp_data
from bot_eo_analyses import repackage_summary_files, \
    run_python_task_or_cl_script
from raft_results_task import raft_results_task

def make_focal_plane_plots():
    #
    # Focal plane heat maps
    #
    fp_configs = {'fe55_BOT_analysis': (('gain', 'psf_sigma'),
                                        ((0.5, 1.5), (3, 5.5)),
                                        ('e-/ADU', 'micron'),
                                        (False, False),
                                        ('1', '1')),
                  'read_noise_BOT': (('read_noise',), ((0, 20),), ('e-',),
                                     (False,), ('1',)),
                  'bright_defects_BOT': (('bright_pixels', 'bright_columns'),
                                         (None, None), (None, None),
                                         (True, True), ('1', '1')),
                  'dark_defects_BOT': (('dark_pixels', 'dark_columns'),
                                       (None, None), (None, None),
                                       (True, True), ('1', '1')),
                  'dark_current_BOT': (('dark_current_95CL',), ((0, 0.3),),
                                       ('e-/pix/s',), (False,), ('1',)),
                  'traps_BOT': (('num_traps',), (None,), (None,), (False,),
                                ('1',)),
                  'cte_BOT': (('cti_high_serial', 'cti_high_parallel',
                               'cti_low_serial', 'cti_low_parallel'),
                              ((0, 1e-5), (0, 1e-5), (0, 1e-5), (0, 1e-5)),
                              (None, None, None, None),
                              (False, False, False, False),
                              ('1e-6', '1e-6', '1e-6', '1e-6')),
                  'flat_pairs_BOT': (('max_frac_dev', 'full_well'),
                                     ((0, 0.05), (0, 2e5)), (None, None),
                                     (False, False), ('1', '1')),
                  'ptc_BOT': (('ptc_gain', 'ptc_a00'), ((0.5, 1.5), (0, 5e-6)),
                              ('e-/ADU', None), (False, False), ('1', '1')),
                  'bright_fatter_BOT': (('bf_xcorr', 'bf_ycorr'), (None, None),
                                        (None, None), (False, False),
                                        ('1', '1'))}

    camera = camera_info.camera_object

    run = siteUtils.getRunNumber()
    unit_id = siteUtils.getUnitId()
    et_results = siteUtils.ETResults(run)
    for schema_name, configs in fp_configs.items():
        for column, z_range, units, use_log10, scale_factor in zip(*configs):
            fig = plt.figure(figsize=(12, 10))
            ax = fig.add_subplot(1, 1, 1)
            try:
                amp_data = et_results.get_amp_data(schema_name, column)
            except KeyError:
                print("No focal plane results for {}: {} in eTraveler."
                      .format(schema_name, column))
            else:
                plot_focal_plane(ax, amp_data, camera=camera, z_range=z_range,
                                 use_log10=use_log10,
                                 scale_factor=scale_factor)
                if units is not None:
                    plt.title('Run {}, {} ({})'.format(run, column, units))
                else:
                    plt.title('Run {}, {}'.format(run, column))
                outfile = '{}_{}_{}.png'.format(unit_id, run, column)
                plt.savefig(outfile)

                # Histogram the amp-level data.
                plt.figure()
                hist_amp_data(amp_data, column, hist_range=z_range,
                              use_log10=use_log10, scale_factor=scale_factor)
                if units is not None:
                    plt.title('Run {}, {} ({})'.format(run, column, units))
                else:
                    plt.title('Run {}, {}'.format(run, column))
                outfile = '{}_{}_{}_hist.png'.format(unit_id, run, column)
                plt.savefig(outfile)

if __name__ == '__main__':
    # Repackage the summary.lims files from all previous harnessed
    # jobs into eotest_results files so that results quantities can be
    # easily retrieved.
    repackage_summary_files()

    raft_results_task_script \
        = os.path.join(os.environ['EOANALYSISJOBSDIR'], 'harnessed_jobs',
                       'raft_results_summary_BOT', 'v0', 'raft_results_task.py')

#    installed_rafts = camera_info.get_installed_raft_names()
    installed_rafts = camera_info.installed_science_rafts
    run_python_task_or_cl_script(raft_results_task, raft_results_task_script,
                                 device_names=installed_rafts)

    make_focal_plane_plots()
