#!/usr/bin/env python
"""
Producer script for raft-level CTE analysis.
"""
from __future__ import print_function
import glob
import lsst.eotest.sensor as sensorTest
import siteUtils
import eotestUtils
from multiprocessor_execution import sensor_analyses

def run_cte_task(sensor_id):
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

    results_file = '%s_eotest_results.fits' % sensor_id
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
                                       % (sensor_id, flux_level)))
        siteUtils.make_png_file(plots.cte_profiles,
                                ('%s_serial_oscan_%s.png' %
                                 (file_prefix, flux_level)),
                                flux_level, sflat_file, mask_files, serial=True)

        siteUtils.make_png_file(plots.cte_profiles,
                                ('%s_parallel_oscan_%s.png' %
                                 (file_prefix, flux_level)),
                                flux_level, sflat_file, mask_files, serial=False)

if __name__ == '__main__':
    sensor_analyses(run_cte_task)
