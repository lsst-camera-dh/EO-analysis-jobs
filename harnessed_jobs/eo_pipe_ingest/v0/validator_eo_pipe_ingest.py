#!/usr/bin/env ipython
"""
Validator script for ingest of eo_pipe results.
"""
import pickle
from lsst.obs.lsst import LsstCam
import lcatr.schema


def get_amp(channel):
    channels = ['C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17',
                'C07', 'C06', 'C05', 'C04', 'C03', 'C02', 'C01', 'C00']
    return channels.index(channel) + 1


schema_names = ('eo_read_noise', 'eo_defects', 'eo_dark_current',
                'eo_eper', 'eo_divisadero_tearing', 'eo_ptc_plots',
                'eo_linearity_plots', 'eo_bf_analysis')
generic_fields = ('schema_name', 'schema_version', 'amp', 'slot', 'raft')

camera = LsstCam.getCamera()

amp_data_file = 'eo_pipe_amp_data.pickle'
with open(amp_data_file, 'rb') as fobj:
    eo_pipe_amp_data = pickle.load(fobj)

results = []
for schema_name in schema_names:
    schema = lcatr.schema.get(schema_name)
    results_fields = [_ for _ in schema if _ not in generic_fields]
    kwargs = {key: None for key in results_fields}
    for det in camera:
        det_name = det.getName()
        kwargs['raft'], kwargs['slot'] = det_name.split('_')
        for amp in det:
            amp_name = amp.getName()
            kwargs['amp'] = get_amp(amp_name)
            for field in results_fields:
                try:
                    kwargs[field] = eo_pipe_amp_data[field][det_name][amp_name]
                except KeyError:
                    pass
            if None in kwargs.values():
                continue
            results.append(lcatr.schema.valid(schema, **kwargs))
lcatr.schema.write_file(results)
lcatr.schema.validate_file()
