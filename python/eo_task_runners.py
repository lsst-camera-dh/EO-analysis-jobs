"""
Functions to run single sensor tasks for raft-level testing.
"""
from __future__ import print_function
import glob
import matplotlib.pyplot as plt
import lsst.eotest.image_utils as imutils
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from correlated_noise import correlated_noise, raft_level_oscan_correlations
import camera_components

__all__ = ['run_fe55_task', 'run_read_noise_task', 'raft_overscan_correlations',
           'raft_overscan_correlations', 'run_bright_pixels_task',
           'run_dark_pixels_task', 'run_trap_task', 'run_dark_current_task',
           'run_cte_task', 'run_flat_pair_task', 'run_ptc_task',
           'run_qe_task']

_default_job_map = {job: siteUtils.getProcessName(job) for job in
                    ('fe55_raft_acq', 'dark_raft_acq', 'sflat_raft_acq',
                     'ppump_raft_acq', 'flat_pair_raft_acq', 'qe_raft_acq')}

def _check_job_map(job_map, *jobs):
    if not job_map:
        return _default_job_map
    missing = []
    for job in jobs:
        if job not in job_map:
            missing.append(job)
        else:
            if job != 'fe55_raft_analysis':
                job_map[job] = siteUtils.getProcessName(job_map[job])
    if missing:
        raise RuntimeError("Expected jobs not found in job_map: %s" %
                           ', '.join(missing))
    return job_map

def getSensorGains(jobname, sensor_id):
    try:
        eotestUtils.getSensorGains(jobname=jobname, sensor_id=sensor_id)
    except RuntimeError:
        data = sensorTest.EOTestResults('%s_eotest_results.fits' % sensor_id)
        return dict(amp, gain for (amp, gain) in zip(data['AMP'], data['GAIN']))

def run_fe55_task(sensor_id, **job_map):
    "Single sensor execution of the Fe55 analysis task."
    job_map = _check_job_map(job_map, 'fe55_raft_acq')
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    fe55_files = siteUtils.dependency_glob('S*/%s_fe55_fe55_*.fits' % sensor_id,
                                           jobname=job_map['fe55_raft_acq'],
                                           description='Fe55 files:')
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                           jobname=job_map['fe55_raft_acq'],
                                           description='Bias files:')
    nf = len(bias_files)
    mean_bias_file = '%(sensor_id)s_mean_bias_%(nf)i.fits' % locals()
    imutils.fits_mean_file(bias_files, mean_bias_file)

    #
    # Create a png zoom of the upper right corner of segment 1 for an Fe55
    # exposure for inclusion in the test report.
    #
    print("processing fe55_zoom:", fe55_files[0])
    siteUtils.make_png_file(sensorTest.fe55_zoom,
                            '%(file_prefix)s_fe55_zoom.png' % locals(),
                            fe55_files[0], size=250, amp=1,
                            annotation='ADU/pixel')

    #
    # Perform analysis of 9-pixel statistics for Fe55 charge clusters.
    #
    try:
        pixel_stats = sensorTest.Fe55PixelStats(fe55_files, sensor_id=sensor_id)

        siteUtils.make_png_file(pixel_stats.pixel_hists,
                                '%s_fe55_p3_p5_hists.png' % file_prefix,
                                pix0='p3', pix1='p5')

        siteUtils.make_png_file(pixel_stats.pixel_diff_profile,
                                '%s_fe55_p3_p5_profiles.png' % file_prefix,
                                pixel_coord='x', pix0='p3', pix1='p5')

        siteUtils.make_png_file(pixel_stats.apflux_profile,
                                '%s_fe55_apflux_serial.png' % file_prefix)

        siteUtils.make_png_file(pixel_stats.apflux_profile,
                                '%s_fe55_apflux_parallel.png' % file_prefix,
                                pixel_coord='y')

    except StandardError as eobj:
        print("Exception raised while creating pixel statistics plots:")
        print(str(eobj))
        print("Skipping these plots.")

    # Roll-off defects mask needs an input file to get the vendor
    # geometry, and will be used for all analyses.
    rolloff_mask_file = '%s_rolloff_defects_mask.fits' % sensor_id
    sensorTest.rolloff_mask(fe55_files[0], rolloff_mask_file)

    task = sensorTest.Fe55Task()
    task.config.temp_set_point = -100.
    task.run(sensor_id, fe55_files, (rolloff_mask_file,), accuracy_req=0.01)

    # Fe55 gain and psf analysis results plots for the test report.
    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    siteUtils.make_png_file(plots.gains,
                            '%s_gains.png' % file_prefix)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_mean_bias.png' % file_prefix,
                            mean_bias_file,
                            title='%s, mean bias frame' % sensor_id,
                            annotation='ADU/pixel, overscan-subtracted')

    fe55_file = glob.glob('%s_psf_results*.fits' % sensor_id)[0]
    siteUtils.make_png_file(plots.fe55_dists,
                            '%s_fe55_dists.png' % file_prefix,
                            fe55_file=fe55_file)

    siteUtils.make_png_file(plots.psf_dists,
                            '%s_psf_dists.png' % file_prefix,
                            fe55_file=fe55_file)


def run_read_noise_task(sensor_id, **job_map):
    """Single sensor execution of the read noise task."""
    job_map = _check_job_map(job_map, 'fe55_raft_acq')
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    bias_files = siteUtils.dependency_glob('S*/%s_fe55_fe55_*.fits' % sensor_id,
                                           jobname=job_map['fe55_raft_acq'],
                                           description='Fe55 files for read noise:')
    gains = getSensorGains('fe55_raft_analysis', sensor_id)

    system_noise = None

    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)

    task = sensorTest.ReadNoiseTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, bias_files, gains, system_noise=system_noise,
             mask_files=mask_files, use_overscan=True)

    # Compute amp-amp correlated noise.
    _, corr_fig, _ = correlated_noise(bias_files, target=0,
                                      make_plots=True, title=sensor_id)
    plt.figure(corr_fig.number)
    plt.savefig('%s_correlated_noise.png' % file_prefix)


def get_bias_files(raft_id=None, **job_map):
    """Get the bias files from the Fe55 acquisition."""
    if not job_map:
        job_map = {'fe55_raft_acq': 'fe55_raft_acq'}
    _check_job_map(job_map, 'fe55_raft_acq')
    if raft_id is None:
        raft_id = os.environ['LCATR_UNIT_ID']
    if not job_map:
        job_map = {'fe55_raft_acq': 'fe55_raft_acq'}
    _check_job_map(job_map, 'fe55_raft_acq')
    raft = camera_components.Raft.create_from_etrav(raft_id)
    bias_files = dict()
    for slot, sensor_id in raft.items():
        bias_files[slot] \
            = siteUtils.dependency_glob('S*/%s_fe55_bias_*.fits' % sensor_id,
                                        jobname=siteUtils.getProcessName('fe55_raft_acq'),
                                        description='Bias files for noise correlations:')[0]
    return bias_files


def raft_overscan_correlations():
    """
    Run the raft_level_oscan_correlations funciton for the raft
    being tested in the current raft-level harnessed job.
    """
    raft_id = os.environ['LCATR_UNIT_ID']
    run = siteUtils.getRunNumber()
    bias_files = get_bias_files(raft_id)
    title = 'Overscan correlations, {}, Run {}'.format(raft_id, run)
    plt.rcParams['figure.figsize'] = (8, 8)
    raft_level_oscan_correlations(bias_files, title=title)
    plt.savefig('{}_{}_overscan_correlations.png'.format(raft_id, run))


def run_bright_pixels_task(sensor_id, **job_map):
    "Single sensor execution of the bright pixels task."
    job_map = _check_job_map(job_map, 'fe55_raft_acq')
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    dark_files = siteUtils.dependency_glob('S*/%s_dark_dark_*.fits' % sensor_id,
                                           jobname=job_map['dark_raft_acq'],
                                           description='Dark files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = getSensorGains('fe55_raft_analysis', sensor_id)

    task = sensorTest.BrightPixelsTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, dark_files, mask_files, gains)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_medianed_dark.png' % file_prefix,
                            '%s_median_dark_bp.fits' % sensor_id,
                            title='%s, medianed dark for bright defects analysis' % sensor_id,
                            annotation='e-/pixel, gain-corrected, bias-subtracted')


def run_dark_pixels_task(sensor_id, **job_map):
    "Single sensor execution of the dark pixels task."
    job_map = _check_job_map(job_map, 'fe55_raft_acq')
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    sflat_files = siteUtils.dependency_glob('S*/%s_sflat_500_flat_H*.fits' % sensor_id,
                                            jobname=job_map['sflat_raft_acq'],
                                            description='Superflat files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)

    task = sensorTest.DarkPixelsTask()
    task.run(sensor_id, sflat_files, mask_files)

    siteUtils.make_png_file(sensorTest.plot_flat,
                            '%s_superflat_dark_defects.png' % file_prefix,
                            '%s_median_sflat.fits' % sensor_id,
                            title='%s, superflat for dark defects analysis' % sensor_id,
                            annotation='ADU/pixel')


def run_trap_task(sensor_id, **job_map):
    job_map = _check_job_map(job_map, 'ppump_raft_acq')
    trap_file = siteUtils.dependency_glob('S*/%s_trap_ppump_*.fits' % sensor_id,
                                          jobname=job_map['ppump_raft_acq'],
                                          description='Trap file:')[0]
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    # Omit rolloff defects mask since a trap in the rolloff edge region can
    # affect the entire column.
    mask_files = [item for item in mask_files
                  if item.find('rolloff_defects') == -1]
    print("Using mask files:")
    for mask_file in mask_files:
        print("  " + mask_file)

    gains = getSensorGains('fe55_raft_analysis', sensor_id)

    task = sensorTest.TrapTask()
    task.run(sensor_id, trap_file, mask_files, gains)


def run_dark_current_task(sensor_id, **job_map):
    "Single sensor execution of dark current analysis."
    job_map = _check_job_map(job_map, 'dark_raft_acq')
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    dark_files = siteUtils.dependency_glob('S*/%s_dark_dark_*.fits' % sensor_id,
                                           jobname=job_map['dark_raft_acq'],
                                           description='Dark files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = getSensorGains('fe55_raft_analysis', sensor_id)

    task = sensorTest.DarkCurrentTask()
    task.config.temp_set_point = -100.
    dark_curr_pixels, dark95s \
        = task.run(sensor_id, dark_files, mask_files, gains)

    results_file \
        = siteUtils.dependency_glob('%s_eotest_results.fits' % sensor_id,
                                    jobname='read_noise_raft')[0]
    eo_results = sensorTest.EOTestResults(results_file)
    read_noise = dict(pair for pair in zip(eo_results['AMP'],
                                           eo_results['TOTAL_NOISE']))

    siteUtils.make_png_file(sensorTest.total_noise_histograms,
                            '%s_total_noise_hists.png' % file_prefix,
                            dark_curr_pixels, read_noise, dark95s,
                            exptime=16, title=sensor_id)

    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)
    siteUtils.make_png_file(plots.total_noise, '%s_noise.png' % file_prefix,
                            dark95s=dark95s)


def run_cte_task(sensor_id):
    "Single sensor execution of the cte task."
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)
    # Omit rolloff defects mask since it would mask some of the edges used
    # in the eper method.
    mask_files = [item for item in mask_files if
                  item.find('rolloff_defects') == -1]
    print("Using mask files:")
    for mask_file in mask_files:
        print("  " + mask_file)

    results_file \
        = siteUtils.dependency_glob('%s_eotest_results.fits' % sensor_id,
                                    jobname='fe55_raft_analysis',
                                    description='Fe55 results file')[0]
    shutil.copy(results_file, os.path.basename(results_file))
    results_file = os.path.basename(results_file)
    sflat_high_files = \
        siteUtils.dependency_glob('S*/%s_sflat_500_flat_H*.fits' % sensor_id,
                                  jobname=siteUtils.getProcessName('sflat_raft_acq'),
                                  description='Superflat high flux files:')

    task = sensorTest.CteTask()
    task.run(sensor_id, sflat_high_files, flux_level='high', gains=gains,
             mask_files=mask_files)

    sflat_low_files = \
        siteUtils.dependency_glob('S*/%s_sflat_500_flat_L*.fits' % sensor_id,
                                  jobname=siteUtils.getProcessName('sflat_raft_acq'),
                                  description='Superflat low flux files:')
    task.run(sensor_id, sflat_low_files, flux_level='low', gains=gains,
             mask_files=mask_files)

    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    superflat_files = sorted(glob.glob('%s_superflat_*.fits' % sensor_id))
    mask_files = [x for x in glob.glob('%s*mask.fits' % sensor_id)
                  if x.find('rolloff') == -1]
    for sflat_file in superflat_files:
        flux_level = 'low'
        if sflat_file.find('high') != -1:
            flux_level = 'high'
        siteUtils.make_png_file(sensorTest.plot_flat,
                                sflat_file.replace('.fits', '.png').replace(sensor_id, file_prefix),
                                sflat_file,
                                title=('%s, CTE supeflat, %s flux '
                                       % (sensor_id, flux_level)),
                                annotation='ADU/pixel')
        siteUtils.make_png_file(plots.cte_profiles,
                                ('%s_serial_oscan_%s.png' %
                                 (file_prefix, flux_level)),
                                flux_level, sflat_file, mask_files, serial=True)

        siteUtils.make_png_file(plots.cte_profiles,
                                ('%s_parallel_oscan_%s.png' %
                                 (file_prefix, flux_level)),
                                flux_level, sflat_file, mask_files, serial=False)


def run_flat_pair_task(sensor_id):
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    flat_files = siteUtils.dependency_glob('S*/%s_flat*flat?_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('flat_pair_raft_acq'),
                                           description='Flat files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    use_exptime = True
    if siteUtils.getSiteName() == 'SLAC':
        # Since slit-width can be set individually for each exposure
        # on TS-8 at IR-2 (LSSTTD-1231), we need to use the MONDIODE
        # keyword for computing the integrated incident flux.
        use_exptime = False

    task = sensorTest.FlatPairTask()
    task.run(sensor_id, flat_files, mask_files, gains,
             linearity_spec_range=(1e4, 9e4), use_exptime=use_exptime)

    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    Ne_bounds = (1e4, 9e4)

    detresp_file = '%s_det_response.fits' % sensor_id
    siteUtils.make_png_file(plots.linearity,
                            '%s_linearity.png' % file_prefix,
                            detresp_file=detresp_file, max_dev=0.03,
                            use_exptime=use_exptime, Ne_bounds=Ne_bounds)
    siteUtils.make_png_file(plots.linearity_resids,
                            '%s_linearity_resids.png' % file_prefix,
                            detresp_file=detresp_file, max_dev=0.03,
                            Ne_bounds=Ne_bounds, use_exptime=use_exptime)


def run_ptc_task(sensor_id):
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    flat_files = siteUtils.dependency_glob('S*/%s_flat*flat?_*.fits' % sensor_id,
                                           jobname=siteUtils.getProcessName('flat_pair_raft_acq'),
                                           description='Flat files:')
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.PtcTask()
    task.run(sensor_id, flat_files, mask_files, gains)

    results_file = '%s_eotest_results.fits' % sensor_id
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)
    siteUtils.make_png_file(plots.ptcs,
                            '%s_ptcs.png' % file_prefix,
                            ptc_file='%s_ptc.fits' % sensor_id)


def run_qe_task(sensor_id):
    "Single sensor execution of the QE task."
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    lambda_files = siteUtils.dependency_glob('S*/%s_lambda_flat_*.fits' % sensor_id,
                                             jobname=siteUtils.getProcessName('qe_raft_acq'),
                                             description='Lambda files:')

    pd_ratio_file = eotestUtils.getPhotodiodeRatioFile()
    if pd_ratio_file is None:
        message = ("The test-stand specific photodiode ratio file is " +
                   "not given in config/%s/eotest_calibrations.cfg."
                   % siteUtils.getSiteName())
        raise RuntimeError(message)

    correction_image = eotestUtils.getIlluminationNonUniformityImage()
    if correction_image is None:
        print()
        print("WARNING: The correction image file is not given in")
        print("config/%s/eotest_calibrations.cfg." % siteUtils.getSiteName())
        print("No correction for non-uniform illumination will be applied.")
        print()
        sys.stdout.flush()

    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.QeTask()
    task.config.temp_set_point = -100.
    task.run(sensor_id, lambda_files, pd_ratio_file, mask_files, gains,
             correction_image=correction_image)

    results_file \
        = siteUtils.dependency_glob('%s_eotest_results.fits' % sensor_id,
                                    jobname='fe55_raft_analysis',
                                    description='Fe55 results file')[0]
    plots = sensorTest.EOTestPlots(sensor_id, results_file=results_file)

    siteUtils.make_png_file(plots.qe,
                            '%s_qe.png' % file_prefix,
                            qe_file='%s_QE.fits' % sensor_id)

    try:
        plots.flat_fields(os.path.dirname(lambda_files[0]),
                          annotation='e-/pixel, gain-corrected, bias-subtracted')
    except StandardError as eobj:
        print("Exception raised while creating flat fields:")
        print(str(eobj))

