#!/usr/bin/env ipython
"""
Validator script for raft-level flat pairs analysis.
"""
import glob
import json
import pickle
import lcatr.schema
import siteUtils
import camera_components
from tearing_detection import persist_tearing_png_files

raft_unit_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_unit_id)

det_map = dict()
results = []
schema = lcatr.schema.get('tearing_detection')
for slot, sensor_id in raft.items():
    det_map[slot] = sensor_id
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

divisidero_plot = glob.glob('*_divisidero.png')[0]
md = dict(DATA_PRODUCT='divisidero_tearing_plot',
          LsstId=raft_unit_id)
results.append(siteUtils.make_fileref(divisidero_plot, metadata=md))

with open(glob.glob('*max_divisidero.json')[0], 'r') as fd:
    max_devs = json.load(fd)

fields = ('max_deviation_10_11',
          'max_deviation_11_12',
          'max_deviation_12_13',
          'max_deviation_13_14',
          'max_deviation_14_15',
          'max_deviation_15_16',
          'max_deviation_16_17',
          'max_deviation_00_01',
          'max_deviation_01_02',
          'max_deviation_02_03',
          'max_deviation_03_04',
          'max_deviation_04_05',
          'max_deviation_05_06',
          'max_deviation_06_07')

divisidero_schema = lcatr.schema.get('divisidero_tearing')
for slot, values in max_devs.items():
    data = {field: max_dev for field, max_dev in zip(fields, values)}
    results.append(lcatr.schema.valid(divisidero_schema, slot=slot,
                                      sensor_id=det_map[slot],
                                      **data))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
