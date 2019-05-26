"""
parsl configuration for running on IR2 diagnostic cluster.
"""
import os
import socket
import logging
import parsl
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider
from parsl.channels import SSHChannel
from parsl.config import Config

__all__ = ['load_ir2_dc_config', 'MAX_PARSL_THREADS']


# Suppress the verbose default debug-level logging:
parsl.set_stream_logger(name='interchange', level=logging.INFO)

SETUP_SCRIPT = os.environ.get('LCATR_SETUP_SCRIPT',
                              os.path.join(os.environ['INST_DIR'], 'setup.sh'))

WORKER_NODE_ADDRESSES = set(['lsst-dc01', 'lsst-dc02', 'lsst-dc03',
                             'lsst-dc04', 'lsst-dc05', 'lsst-dc06',
                             'lsst-dc07', 'lsst-dc08', 'lsst-dc09',
                             'lsst-dc10'])
MOTHER_NODE_ADDRESS = socket.gethostname().split('.')[0]
WORKER_NODE_ADDRESSES.remove(MOTHER_NODE_ADDRESS)

NCORES = 28
MAX_PARSL_THREADS = len(WORKER_NODE_ADDRESSES) * NCORES


def script_dir(hostname, root_dir='.'):
    """
    Compose a directory path for local provider scripts.
    """
    return os.path.join(os.path.abspath(root_dir), 'runinfo', hostname)


def load_ir2_dc_config():
    """
    Load the parsl config for ad-hoc providers.
    """
    try:
        parsl.DataFlowKernelLoader.dfk()
        print("parsl config is already loaded.")
        return
    except RuntimeError:
        pass

    executors = []

    for host in WORKER_NODE_ADDRESSES:
        channel = SSHChannel(hostname=host, script_dir=script_dir(host))
        provider = LocalProvider(channel=channel,
                                 worker_init='source %s' % SETUP_SCRIPT)
        executors.append(HighThroughputExecutor(label=host,
                                                address=MOTHER_NODE_ADDRESS,
                                                worker_debug=False,
                                                provider=provider))

    config = Config(executors=executors, strategy=None)

    parsl.load(config)
