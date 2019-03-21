#!/usr/bin/env python
"""
Validator script for BOT bias frame generation.
"""
import glob
import lcatr.schema
import siteUtils
from camera_components import camera_info
from bot_eo_analyses import make_file_prefix

det_names = camera_info.get_det_names()

run = siteUtils.getRunNumber()

results = []
for det_name in det_names:
    file_prefix = make_file_prefix(run, det_name)
    bias_frame = glob.glob('{}_median_bias.fits'.format(file_prefix))[0]
    results.append(siteUtils.make_fileref(bias_frame))

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
