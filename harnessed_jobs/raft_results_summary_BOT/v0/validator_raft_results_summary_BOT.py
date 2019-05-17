#!/usr/bin/env python
"""
Validator script for BOT raft-level results summaries.
"""
import glob
import lcatr.schema
import siteUtils
from camera_components import camera_info
from bot_eo_validators import validate_raft_results

results = []
results = validate_raft_results(results, camera_info.get_raft_names())

#
# Validate focal plane heat map plots
#
unit_id = siteUtils.getUnitId()
run = siteUtils.getRunNumber()
png_files = glob.glob('{}_{}*.png'.format(unit_id, run))
results.extend(siteUtils.persist_png_files('', unit_id, png_files=png_files))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
