#!/usr/bin/env python
"""
Validator script for raft-level CTE analysis.
"""
import glob
import lsst.eotest.sensor as sensorTest
import lcatr.schema
import siteUtils
import eotestUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    superflats = glob.glob('%(sensor_id)s_superflat_*.fits' % locals())
    for item in superflats:
        eotestUtils.addHeaderData(item, FILENAME=item,
                                  DATE=eotestUtils.utc_now_isoformat())
    results = [siteUtils.make_fileref(x, folder=slot) for x in superflats]

    results_file = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(results_file)
    amps = data['AMP']

    cti_high_serial = data['CTI_HIGH_SERIAL']
    cti_high_serial_error = data['CTI_HIGH_SERIAL_ERROR']
    cti_high_parallel = data['CTI_HIGH_PARALLEL']
    cti_high_parallel_error = data['CTI_HIGH_PARALLEL_ERROR']

    cti_low_serial = data['CTI_LOW_SERIAL']
    cti_low_serial_error = data['CTI_LOW_SERIAL_ERROR']
    cti_low_parallel = data['CTI_LOW_PARALLEL']
    cti_low_parallel_error = data['CTI_LOW_PARALLEL_ERROR']

    for values in zip(amps,
                      cti_high_serial, cti_high_serial_error,
                      cti_high_parallel, cti_high_parallel_error,
                      cti_low_serial, cti_low_serial_error,
                      cti_low_parallel, cti_low_parallel_error):
        results.append(lcatr.schema.valid(lcatr.schema.get('cte'),
                                          amp=values[0],
                                          cti_high_serial=values[1],
                                          cti_high_serial_error=values[2],
                                          cti_high_parallel=values[3],
                                          cti_high_parallel_error=values[4],
                                          cti_low_serial=values[5],
                                          cti_low_serial_error=values[6],
                                          cti_low_parallel=values[7],
                                          cti_low_parallel_error=values[8],
                                          slot=slot,
                                          sensor_id=sensor_id))

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
