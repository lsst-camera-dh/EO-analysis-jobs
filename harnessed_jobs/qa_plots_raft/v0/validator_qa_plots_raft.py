#!/usr/bin/env python
import lcatr.schema
import siteUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()
    results.extend(siteUtils.persist_png_files('%s*.png' % sensor_id,
                                               ccd_vendor, sensor_id,
                                               'QA', 'EO', folder=slot))
results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
