#!/usr/bin/env ipython
"""
Producer script for BOT raft-level results summaries.
"""
def raft_results_task(raft_name):
    """Task to aggregate data for raft-level plots and results."""
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import lsst.eotest.sensor as sensorTest
    import lsst.eotest.raft as raftTest
    import siteUtils
    from camera_components import camera_info
    from bot_eo_analyses import get_raft_files_by_slot, make_file_prefix,\
        get_amplifier_gains

    def plt_savefig(filename):
        plt.savefig(filename)
        plt.close()

    # Get results files for each CCD in the raft.
    try:
        results_files \
            = get_raft_files_by_slot(raft_name, 'eotest_results.fits')
    except FileNotFoundError:
        print("No raft-level results for", raft_name)
        return None

    # Determine the total number of pixels and number of edge rolloff
    # pixels for the types of CCDs in this raft and update the results
    # files.  This info will be used in computing the pixel defect
    # compliance.  Use one of the median bias files for this since they
    # should be available no matter which analysis tasks are run.
    bias_frames = get_raft_files_by_slot(raft_name, 'median_bias.fits',
                                         jobname='bias_frame_BOT')
    try:
        mask_files = get_raft_files_by_slot(raft_name,
                                            'edge_rolloff_mask.fits')
    except FileNotFoundError:
        input_mask = None
    else:
        input_mask = list(mask_files.values())[0]
    total_num, rolloff_mask \
        = sensorTest.pixel_counts(list(bias_frames.values())[0],
                                  input_mask=input_mask)

    # Exposure time (in seconds) for 95th percentile dark current shot
    # noise calculation.
    exptime = 15.

    # Update the eotest results files.
    for filename in results_files.values():
        eotest_results = sensorTest.EOTestResults(filename)
        eotest_results.add_ccd_result('TOTAL_NUM_PIXELS', total_num)
        eotest_results.add_ccd_result('ROLLOFF_MASK_PIXELS', rolloff_mask)
        shot_noise = eotest_results['DARK_CURRENT_95']*exptime
        total_noise = np.sqrt(eotest_results['READ_NOISE']**2 + shot_noise)
        add_max_frac_dev = 'MAX_FRAC_DEV' not in eotest_results.colnames
        for i, amp in enumerate(eotest_results['AMP']):
            if add_max_frac_dev:
                eotest_results.add_seg_result(amp, 'MAX_FRAC_DEV', 0.)
            eotest_results.add_seg_result(amp, 'DC95_SHOT_NOISE',
                                          np.float(shot_noise[i]))
            try:
                eotest_results['TOTAL_NOISE'][i] = total_noise[i]
            except KeyError:
                eotest_results.add_seg_result(amp, 'TOTAL_NOISE',
                                              np.float(total_noise[i]))
        eotest_results.write(filename)

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, raft_name)
    title = '{}, {}'.format(run, raft_name)

    gains = {slot_name: get_amplifier_gains(results_files[slot_name])
             for slot_name in results_files}

    # Update the gains in the results files with the retrieved values.
    for slot_name, ccd_gains in gains.items():
        results = sensorTest.EOTestResults(results_files[slot_name])
        for amp, gain in ccd_gains.items():
            results.add_seg_result(amp, 'GAIN', gain)
        results.write()

    # Extract dark currents for each amplifier in the raft.
    dark_currents = dict()
    for slot_name, results_file in results_files.items():
        results = sensorTest.EOTestResults(results_file)
        try:
            dark_currents[slot_name] \
                = dict(_ for _ in zip(results['AMP'],
                                      results['DARK_CURRENT_MEDIAN']))
        except KeyError:
            dark_currents[slot_name] = dict({amp: 0 for amp in range(1, 17)})

    png_files = []
    # Median bias mosaic
    median_bias = raftTest.make_raft_mosaic(bias_frames, bias_subtract=False)
    median_bias.plot(title='%s, median bias frames' % title,
                     annotation='ADU/pixel', rotate=180)
    png_files.append('{}_median_bias.png'.format(file_prefix))
    plt_savefig(png_files[-1])
    del median_bias

    # Dark mosaic
    try:
        dark_files = get_raft_files_by_slot(raft_name, 'median_dark_bp.fits')
    except FileNotFoundError as eobj:
        print(eobj)
    else:
        dark_mosaic = raftTest.make_raft_mosaic(dark_files, gains=gains,
                                                bias_frames=bias_frames)
        dark_mosaic.plot(title='{}, medianed dark frames'.format(title),
                         annotation='e-/pixel, gain-corrected, bias-subtracted',
                         rotate=180)
        png_files.append('{}_medianed_dark.png'.format(file_prefix))
        plt_savefig(png_files[-1])
        del dark_mosaic

    # High flux superflat mosaic.
    try:
        sflat_high_files \
            = get_raft_files_by_slot(raft_name, 'superflat_high.fits')
    except FileNotFoundError as eobj:
        print(eobj)
    else:
        sflat_high = raftTest.make_raft_mosaic(sflat_high_files, gains=gains,
                                               bias_frames=bias_frames,
                                               dark_currents=dark_currents)
        sflat_high.plot(title='%s, high flux superflat' % title,
                        annotation='e-/pixel, gain-corrected, bias-subtracted',
                        rotate=180)
        png_files.append('{}_superflat_high.png'.format(file_prefix))
        plt_savefig(png_files[-1])
        del sflat_high

    # Low flux superflat mosaic.
    try:
        sflat_low_files \
            = get_raft_files_by_slot(raft_name, 'superflat_low.fits')
    except FileNotFoundError as eobj:
        print(eobj)
    else:
        sflat_low = raftTest.make_raft_mosaic(sflat_low_files, gains=gains,
                                              bias_frames=bias_frames,
                                              dark_currents=dark_currents)
        sflat_low.plot(title='%s, low flux superflat' % title,
                       annotation='e-/pixel, gain-corrected, bias-subtracted',
                       rotate=180)
        png_files.append('{}_superflat_low.png'.format(file_prefix))
        plt_savefig(png_files[-1])
        del sflat_low

    # QE images at various wavelengths and filters
    acq_jobname = siteUtils.getProcessName('BOT_acq')
    for wl in ('SDSSu', 'SDSSg', 'SDSSr', 'SDSSi', 'SDSSz', 'SDSSY',
               '480nm', '650nm', '750nm', '870nm', '950nm', '970nm'):
        print("Processing %s image" % wl)
        pattern = 'lambda_flat_{}*/*_{}_*.fits'.format(wl, raft_name)
        print(pattern)
        print(acq_jobname)
        files = siteUtils.dependency_glob(pattern, acq_jobname=acq_jobname)
        if not files:
            print("no files found")
            continue
        lambda_files = dict()
        for item in files:
            slot_name = os.path.basename(item).split('_')[-1].split('.')[0]
            lambda_files[slot_name] = item
        flat = raftTest.make_raft_mosaic(lambda_files, gains=gains,
                                         bias_frames=bias_frames,
                                         dark_currents=dark_currents)
        flat.plot(title='%s, %s' % (title, wl),
                  annotation='e-/pixel, gain-corrected, bias-subtracted',
                  rotate=180)
        png_files.append('{}_{}_flat.png'.format(file_prefix, wl))
        plt_savefig(png_files[-1])
        del flat

    # TODO: QE summary plot

    # Plots of read noise, nonlinearity, serial and parallel CTI,
    # PSF size, and gains from Fe55 and PTC.
    spec_plots = raftTest.RaftSpecPlots(results_files)
    columns = 'READ_NOISE DC95_SHOT_NOISE TOTAL_NOISE'.split()
    try:
        spec_plots.make_multi_column_plot(columns, 'noise per pixel (-e rms)',
                                          spec=9, title=title,
                                          ybounds=(-1, 100))
        png_files.append('%s_total_noise.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        pass

    try:
        spec_plots.make_plot('MAX_FRAC_DEV',
                             'non-linearity (max. fractional deviation)',
                             spec=0.03, title=title, ybounds=(0, 0.1))
        png_files.append('%s_linearity.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        pass

    try:
        spec_plots.make_multi_column_plot(('CTI_LOW_SERIAL', 'CTI_HIGH_SERIAL'),
                                          'Serial CTI (ppm)', spec=(5e-6, 3e-5),
                                          title=title, yscaling=1e6, yerrors=True,
                                          colors='br', ybounds=(-1e-5, 6e-5))
        png_files.append('%s_serial_cti.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        pass

    try:
        spec_plots.make_multi_column_plot(('CTI_LOW_PARALLEL', 'CTI_HIGH_PARALLEL'),
                                          'Parallel CTI (ppm)', spec=3e-6,
                                          title=title, yscaling=1e6, yerrors=True,
                                          colors='br', ybounds=(-1e-5, 6e-5))
        png_files.append('%s_parallel_cti.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        pass

    try:
        spec_plots.make_plot('PSF_SIGMA', 'PSF sigma (microns)', spec=5.,
                             title=title, ybounds=(0, 5.2))
        png_files.append('%s_psf_sigma.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        # PSF_SIGMA not available so skip this plot
        pass

    try:
        spec_plots.make_multi_column_plot(('GAIN', 'PTC_GAIN'),
                                          'System Gain (e-/ADU)',
                                          yerrors=True, title=title,
                                          colors='br', ybounds=(0, 3))
        png_files.append('%s_system_gain.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        # PTC_GAIN data not available so skip this plot.
        pass

    try:
        spec_plots.make_plot('DARK_CURRENT_95',
                             '95th percentile dark current (e-/pixel/s)',
                             spec=0.2, title=title, ybounds=(-0.01, 1))
        png_files.append('%s_dark_current.png' % file_prefix)
        plt_savefig(png_files[-1])
    except KeyError:
        pass

    # Make bias frame stats time history plots for the current raft.
    pattern = f'{raft_name}_{run}_bias_frame_stats.pkl'
    try:
        stats_file = siteUtils.dependency_glob(pattern,
                                               jobname='bias_frame_BOT')[0]
    except IndexError:
        pass
    else:
        file_prefix = make_file_prefix(run, raft_name)
        df_raft = pd.read_pickle(stats_file)
        slots = '20 21 22 10 11 12 00 01 02'.split()

        fig = plt.figure(figsize=(16, 16))
        for i, slot in enumerate(slots, 1):
            fig.add_subplot(3, 3, i)
            df = df_raft.query(f'slot == "{slot}"')
            amps = sorted(list(set(df['amp'])))
            for amp in amps:
                my_df = df.query(f'amp == {amp}')
                plt.scatter(my_df['tseqnum'], my_df['mean'], s=2, label=f'amp')
            xmin, xmax, _, _ = plt.axis()
            plt.xlim(xmin, 1.2*(xmax - xmin) + xmin)
            plt.legend(fontsize='x-small')
            plt.xlabel('test sequence number')
            plt.ylabel('mean signal (ADU)')
            plt.title(slot)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(f'{file_prefix}, bias stability')
        png_file = f'{file_prefix}_bias_stability_mean.png'
        png_files.append(png_file)
        plt_savefig(png_file)

        fig = plt.figure(figsize=(16, 16))
        for i, slot in enumerate(slots, 1):
            fig.add_subplot(3, 3, i)
            df = df_raft.query(f'slot == "{slot}"')
            amps = sorted(list(set(df['amp'])))
            for amp in amps:
                my_df = df.query(f'amp == {amp}')
                plt.scatter(my_df['tseqnum'], my_df['stdev'], s=2, label=f'amp')
            xmin, xmax, _, _ = plt.axis()
            plt.xlim(xmin, 1.2*(xmax - xmin) + xmin)
            plt.legend(fontsize='x-small')
            plt.xlabel('test sequence number')
            plt.ylabel('stdev (ADU)')
            plt.title(slot)
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(f'{file_prefix}, bias stability')
        png_file = f'{file_prefix}_bias_stability_stdev.png'
        png_files.append(png_file)
        plt_savefig(png_file)

    png_file_list = '{}_raft_results_task_png_files.txt'.format(raft_name)
    with open(png_file_list, 'w') as output:
        for item in png_files:
            if os.path.isfile(item):
                output.write('{}\n'.format(item))

    return None

if __name__ == '__main__':
    import sys
    raft_name = sys.argv[1]
    raft_results_task(raft_name)
