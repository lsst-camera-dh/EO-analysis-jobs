
import os

import parsl
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider
from parsl.channels import SSHChannel
from parsl.config import Config


SETUP_SCRIPT = '/home/echarles/lsst_setup.bashrc'

WORKER_NODE_ADDRESSES = ['lsst-dc01', 'lsst-dc02', 'lsst-dc03',
                         'lsst-dc04', 'lsst-dc05', 'lsst-dc06',
                         'lsst-dc07', 'lsst-dc08', 'lsst-dc09']

MOTHER_NODE_ADDRESS = 'lsst-dc10'
NCORES = 48

MAX_PARSL_THREADS = len(WORKER_NODE_ADDRESSES) * NCORES

PARSL_LOADED = False

def load_ir2_dc_config():

    providers = [LocalProvider(channel=SSHChannel(hostname=hostn,
                                                  no_auth=True),
                               worker_init='source %s' % SETUP_SCRIPT) for hostn in WORKER_NODE_ADDRESSES]

    executors = [HighThroughputExecutor(label="Ad-Hoc_%i" % iprov,
                                        address=MOTHER_NODE_ADDRESS,
                                        worker_debug=True,
                                        cores_per_worker=NCORES,
                                        provider=provider) for iprov, provider in enumerate(providers)]

    config = Config(executors=executors, strategy=None)
    
    parsl.load(config)

    PARSL_LOADED = True
