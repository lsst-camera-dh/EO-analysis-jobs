#!/usr/bin/env python
"""
Validator script for raft-level QE analysis.
"""
import glob
from collections import OrderedDict
import astropy.io.fits as fits
import numpy as np
import lcatr.schema
import siteUtils
import eotestUtils
import camera_components

raft_id = siteUtils.getUnitId()
raft = camera_components.Raft.create_from_etrav(raft_id)

results = []
for slot, sensor_id in raft.items():
    ccd_vendor = sensor_id.split('-')[0].upper()

    qe_data = fits.open('%s_QE.fits' % sensor_id)['QE_BANDS'].data
    QE = OrderedDict((band, []) for band in qe_data.field('BAND'))
    for amp in range(1, 17):
        values = qe_data.field('AMP%02i' % amp)
        for band, value in zip(QE, values):
            QE[band].append(value)

    for band in QE:
        results.append(lcatr.schema.valid(lcatr.schema.get('qe_raft_analysis'),
                                          band=band, QE=np.mean(QE[band]),
                                          slot=slot, sensor_id=sensor_id))

    qe_acq_job_id = siteUtils.get_prerequisite_job_id('S*/%s_lambda_flat_*.fits' % sensor_id,
                                                      jobname='qe_raft_acq_sim')
    md = dict(photodiode_ratio_file=dict(JOB_ID=qe_acq_job_id),
              illumination_non_uniformity_file=dict(JOB_ID=qe_acq_job_id))
    results.extend(eotestUtils.eotestCalibsPersist('photodiode_ratio_file',
                                                   'illumination_non_uniformity_file',
                                                   metadata=md))
    qe_files = glob.glob('*QE*.*')
    for item in qe_files:
        if item.endswith('.fits'):
            eotestUtils.addHeaderData(item, LSST_NUM=sensor_id,
                                      TESTTYPE='LAMBDA',
                                      DATE=eotestUtils.utc_now_isoformat(),
                                      CCD_MANU=ccd_vendor)
    results.extend([siteUtils.make_fileref(item, folder=slot)
                    for item in qe_files])

results.extend(siteUtils.jobInfo())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
