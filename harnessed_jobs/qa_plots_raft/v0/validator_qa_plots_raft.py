#!/usr/bin/env python
"""
Validator script for qa_plots.
"""
import lcatr.schema
import siteUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

metadata = dict(TESTTYPE='QA', TEST_CATEGORY='EO')
results = siteUtils.persist_png_files('%s*.png' % raft_id,
                                      raft_id, metadata=metadata)

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
