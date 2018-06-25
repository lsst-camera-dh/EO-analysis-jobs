#!/usr/bin/env python
"""
Validator script for collect_raft_results job.
"""
from __future__ import print_function
import pickle
import lcatr.schema
import siteUtils
import eotestUtils
import camera_components

results = []

run_number = siteUtils.getRunNumber()

md = siteUtils.DataCatalogMetadata(ORIGIN=siteUtils.getSiteName(),
                                   TEST_CATEGORY='EO',
                                   DATA_PRODUCT='EOTEST_RESULTS')

# Persist eotest_results files for each sensor.
raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()
    results_file = '%s_eotest_results.fits' % sensor_id
    eotestUtils.addHeaderData(results_file, LSST_NUM=sensor_id,
                              DATE=eotestUtils.utc_now_isoformat(),
                              CCD_MANU=ccd_vendor,
                              RUNNUM=run_number)
    results.append(lcatr.schema.fileref.make(results_file,
                                             metadata=md(CCD_MANU=ccd_vendor,
                                                         LSST_NUM=sensor_id,
                                                         SLOT=slot,
                                                         LsstId=raft_id)))

# Persist the png files.
metadata = dict(CCD_MANU=ccd_vendor, TEST_CATEGORY='EO')
results.extend(siteUtils.persist_png_files('%s*.png' % raft_id,
                                           raft_id, metadata=metadata))

# Persist the tearing results.
file_prefix = '%s_%s' % (raft_id, run_number)
with open('%s_tearing_stats.pkl' % file_prefix, 'rb') as input_:
    tearing_stats = pickle.load(input_)
schema = lcatr.schema.get('tearing_detection')
for values in tearing_stats:
    stats = dict(kv_pair for kv_pair in
                 zip(('job_name', 'subset', 'slot', 'sensor_id', 'detections'),
                     values))
    results.append(lcatr.schema.valid(schema, **stats))

results.extend(siteUtils.jobInfo())
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
