#!/usr/bin/env ipython
"""
Producer script for raft-level traps analysis.
"""

def run_trap_task(sensor_id):
    """Run the traps analysis for a single sensor."""
    import lsst.eotest.sensor as sensorTest
    import siteUtils
    import eotestUtils

    trap_file = siteUtils.dependency_glob('S*/%s_trap_ppump_*.fits' % sensor_id,
                                          jobname=siteUtils.getProcessName('ppump_raft_acq'),
                                          description='Trap file:')[0]
    mask_files = \
        eotestUtils.glob_mask_files(pattern='%s_*mask.fits' % sensor_id)
    # Omit rolloff defects mask since a trap in the rolloff edge region can
    # affect the entire column.
    bias_frame = siteUtils.dependency_glob('%s_sflat*median_bias.fits'
                                           % sensor_id,
                                           description='Super bias frame:')[0]
    mask_files = [item for item in mask_files
                  if item.find('rolloff_defects') == -1]
    print("Using mask files:")
    for mask_file in mask_files:
        print("  " + mask_file)

    gains = eotestUtils.getSensorGains(jobname='fe55_raft_analysis',
                                       sensor_id=sensor_id)

    task = sensorTest.TrapTask()
    task.run(sensor_id, trap_file, mask_files, gains, bias_frame=bias_frame)

if __name__ == '__main__':
    import os
    from multiprocessor_execution import sensor_analyses

    processes = 9                # Reserve 1 process per CCD.
    if os.environ.get('LCATR_SKIP_TRAPS_ANALYSIS', 'False') == 'True':
        pass
    else:
        sensor_analyses(run_trap_task, processes=processes)
