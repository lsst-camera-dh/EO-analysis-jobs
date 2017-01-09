#!/usr/bin/env python
"""
Validator script for raft-level PRNU analysis.
"""
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
    results_file = '%s_eotest_results.fits' % sensor_id
    prnu_results = fits.open(results_file)['PRNU_RESULTS'].data

    for wl, stdev, mean in zip(prnu_results['WAVELENGTH'],
                               prnu_results['STDEV'], prnu_results['MEAN']):
        results.append(lcatr.schema.valid(lcatr.schema.get('prnu'),
                                          wavelength=int(np.round(wl)),
                                          pixel_stdev=stdev, pixel_mean=mean,
                                          slot=slot,
                                          sensor_id=sensor_id))

        qe_acq_job_id = \
            siteUtils.get_prerequisite_job_id('S*/%s_lambda_flat_*.fits' % sensor_id,
                                              jobname='qe_raft_acq_sim')
        md = dict(illumination_non_uniformity_file=dict(JOB_ID=qe_acq_job_id))
        results.extend(eotestUtils.eotestCalibsPersist('illumination_non_uniformity_file',
                                                       metadata=md))

results.extend(siteUtils.jobInfo())
results.append(eotestUtils.eotestCalibrations())

lcatr.schema.write_file(results)
lcatr.schema.validate_file()
