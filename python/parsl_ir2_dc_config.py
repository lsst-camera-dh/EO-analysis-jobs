
import parsl
import os

from parsl.app.app import python_app
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider
from parsl.channels import SSHChannel
from parsl.config import Config


SETUP_SCRIPT = '/home/echarles/lsst_setup.bashrc'

#WORKER_NODE_ADDRESSES = ['lsst-dc01', 'lsst-dc02', 'lsst-dc03',
#                         'lsst-dc04', 'lsst-dc05', 'lsst-dc06',
#                         'lsst-dc07', 'lsst-dc08', 'lsst-dc09']

WORKER_NODE_ADDRESSES = ['lsst-dc02', 'lsst-dc03']

MOTHER_NODE_ADDRESS = 'lsst-dc10'
NCORES = 48

PROVIDERS = [LocalProvider(channel=SSHChannel(hostname=hostn,
                                              no_auth=True),
                           worker_init='source %s' % SETUP_SCRIPT) for hostn in WORKER_NODE_ADDRESSES]

EXECUTORS = [HighThroughputExecutor(label="Ad-Hoc_%i" % iprov,
                                    address=MOTHER_NODE_ADDRESS,
                                    worker_debug=True,
                                    cores_per_worker=NCORES,
                                    provider=provider) for iprov, provider in enumerate(PROVIDERS)]

IR2_DC_CONFIG = Config(executors=EXECUTORS, strategy=None)

MAX_PARSL_THREADS = len(WORKER_NODE_ADDRESSES) * NCORES
