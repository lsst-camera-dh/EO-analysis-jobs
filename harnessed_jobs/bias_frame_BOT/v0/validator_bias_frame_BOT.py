#!/usr/bin/env python
"""
Validator script for BOT bias frame generation.
"""
import lcatr.schema
import siteUtils
from camera_components import camera_info
from bot_eo_validators import validate_bias_frames

results = []
results = validate_bias_frames(results, camera_info.get_det_names())
results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
