#!/usr/bin/env python
"""
Validator script for BOT raft-level results summaries.
"""
import lcatr.schema
import siteUtils
from camera_components import camera_info
from bot_eo_validators import validate_raft_results

results = []
results = validate_raft_results(results, camera_info.get_raft_names())
results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
