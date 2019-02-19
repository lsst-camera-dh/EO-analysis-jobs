#!/usr/bin/env python
"""
Validator script for raft-level flat pairs analysis.
"""
import glob
import pickle
import lcatr.schema
import siteUtils
import camera_components
from tearing_detection import persist_tearing_png_files

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
schema = lcatr.schema.get('tearing_detection')
for slot, sensor_id in raft.items():
    file_prefix = '%s_%s' % (sensor_id, siteUtils.getRunNumber())
    with open('%s_tearing_stats.pkl' % file_prefix, 'rb') as input_:
        tearing_stats = pickle.load(input_)
    for values in tearing_stats:
        stats = dict(_ for _ in zip(('job_name', 'subset', 'sensor_id',
                                     'detections', 'slot'),
                                    list(values) + [slot]))
        results.append(lcatr.schema.valid(schema, **stats))

png_files = sorted(glob.glob('*_tearing.png'))
results.extend(persist_tearing_png_files(png_files))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
