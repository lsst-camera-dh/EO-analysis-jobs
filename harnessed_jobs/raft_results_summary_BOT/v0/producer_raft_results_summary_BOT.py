#!/usr/bin/env ipython
"""
Producer script for BOT raft-level results summaries.
"""
import os
import subprocess
import matplotlib.pyplot as plt
import siteUtils
from camera_components import camera_info
from focal_plane_plotting import plot_focal_plane, hist_amp_data
from bot_eo_analyses import repackage_summary_files, \
    run_python_task_or_cl_script
from raft_results_task import raft_results_task
from stage_bot_data import clean_up_scratch

def make_focal_plane_plots():
    #
    # Focal plane heat maps
    #
    fp_configs = {'fe55_BOT_analysis': (('gain', 'psf_sigma'),
                                        ((0.5, 1.6), (3, 5.5)),
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
                  'flat_pairs_BOT': (('max_frac_dev', 'max_observed_signal',
                                      'row_mean_var_slope',
                                      'linearity_turnoff'),
                                     ((0, 0.05), (5e4, 2.e5), (0, 2),
                                      (5e4, 2.e5)),
                                     (None, 'ADU', None, 'ADU'),
                                     (False, False, False, False),
                                     ('1', '1', '1', '1')),
                  'ptc_BOT': (('ptc_gain', 'ptc_a00', 'ptc_turnoff'),
                              ((0.5, 1.6), (0, 5e-6), (5e4, 2e5)),
                              ('e-/ADU', None, 'ADU'),
                              (False, False, False), ('1', '1', '1')),
                  'brighter_fatter_BOT': (('bf_xcorr', 'bf_ycorr'),
                                          ((-5, 5), (-5, 5)),
                                          (None, None), (False, False),
                                          ('1', '1')),
                  'tearing_stats_BOT': (('tearing_detections',), (None,),
                                        (None,), (False,), ('1',)),
                  'divisadero_tearing_BOT': (('divisadero_max_dev',),
                                             ((0, 0.05),), (None,), (False,),
                                             ('1',))}

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

    installed_rafts = camera_info.get_installed_raft_names()
    run_python_task_or_cl_script(raft_results_task, raft_results_task_script,
                                 device_names=installed_rafts)

    make_focal_plane_plots()

    # Make BOT EO report static pages.
    jh_stage_dir = os.path.join(os.environ['LCATR_STAGE_ROOT'],
                                os.environ['LCATR_UNIT_TYPE'],
                                os.environ['LCATR_UNIT_ID'])
    run_number = os.environ.get('LCATR_RUN_NUMBER', '0000')
    run_dir = os.path.join(jh_stage_dir, run_number)

    # Set-up symlinked directories to accommodate eo_task.py assumptions
    # about folder locations and names.
    indir = 'indir'
    bot_folder = 'bot'
    os.makedirs(indir, exist_ok=True)
    dest_dir = os.path.join(indir, bot_folder)
    if not os.path.isdir(run_dir):
        # This traveler is being run using the fake_eT.py script, so
        # we need to modify the destination directory for the symlink
        # to include the run number under an actual `bot_folder`
        # sub-directory.
        os.makedirs(dest_dir)
        dest_dir = os.path.join(dest_dir, run_number)
    os.symlink(jh_stage_dir, dest_dir)

    # Locations of template files for the static reports.
    report_template = os.path.join(os.environ['EO_UTILITIES_DIR'], 'templates',
                                   'eo_html_report.yaml')
    css_file = os.path.join(os.environ['EO_UTILITIES_DIR'], 'templates',
                            'style.css')

    # Write static html to the BOT_EO_reports area, using the fs3
    # location by default.
    htmldir = os.environ.get('LCATR_HTML_DIR',
                             '/gpfs/slac/lsst/fs3/g/data/BOT_EO_reports')
    command = (f'eo_task.py ReportRun --template_file {report_template} '
               f'--indir {indir} --plot_report_action copy '
               f'--runs {run_number} '
               f'--css_file {css_file} --htmldir {htmldir} --overwrite')
    print(command)
    subprocess.check_call(command, shell=True)

    # Clean up scratch areas
    clean_up_scratch(run_number)
