#!/usr/bin/env python
"""
Validator script for raft-level dark current analysis.
"""
import lsst.eotest.sensor as sensorTest
import lcatr.schema
import siteUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    results_file = '%s_eotest_results.fits' % sensor_id
    data = sensorTest.EOTestResults(results_file)

    amps = data['AMP']
    dc95s = data['DARK_CURRENT_95']
    for amp, dc95 in zip(amps, dc95s):
        results.append(lcatr.schema.valid(lcatr.schema.get('dark_current_raft'),
                                          amp=amp, dark_current_95CL=dc95,
                                          slot=slot,
                                          sensor_id=sensor_id))

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
