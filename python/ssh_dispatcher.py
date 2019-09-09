"""
Use ssh to run task scripts on remote nodes that are connected to
a common file system and keep track of remote tasks via log files.
"""
import os
import glob
import copy
import time
import socket
import logging
import subprocess
import multiprocessing
import numpy as np

__all__ = ['ssh_device_analysis_pool']

logging.basicConfig(format='%(asctime)s %(name)s: %(message)s')

class TaskRunner:
    def __init__(self, script, working_dir, setup, verbose=False):
        self.params = script, working_dir, setup
        self.verbose = verbose
        self.log_dir = os.path.join(working_dir, 'logging')
        os.makedirs(self.log_dir, exist_ok=True)

    def log_file(self, task_id, clean_up=True):
        """
        Create a log filename from the task_id and clean up any
        existing log files in the working directory.
        """
        script, working_dir, _ = self.params
        task_name = os.path.basename(script).split('.')[0]
        my_log_file = os.path.join(self.log_dir, f'{task_name}_{task_id}.log')
        if clean_up and os.path.isfile(my_log_file):
            os.remove(my_log_file)
        return my_log_file

    def __call__(self, remote_host, log_file, *args):
        """
        Function call-back for launching the remote process via ssh.
        """
        logger = logging.getLogger('TaskRunner.__call__')
        logger.setLevel(logging.INFO)

        script, working_dir, setup = self.params
        task_id = args[0]
        command = f'ssh {remote_host} '
        command += f'"cd {working_dir}; source {setup}; ({script} '
        command += ' '.join([str(_) for _ in args])
        command += ' && echo Task succeeded)'
        command += f' >& {log_file}&"'
        if self.verbose:
            logger.info(command)
        subprocess.check_call(command, shell=True)
        return os.path.join(working_dir, log_file)

    @staticmethod
    def monitor_tasks(log_files, max_time=None, interval=1, logger=None):
        """
        Static function to keep track of remote processes by
        inspecting log files.
        """
        if logger is None:
            logger = logging.getLogger('TaskRunner.monitor_tasks')
            logger.setLevel(logging.INFO)

        my_log_files = copy.deepcopy(log_files)
        # Poll log files for completion of each task:
        t0 = time.time()
        while my_log_files:
            if max_time is not None and time.time() - t0 > max_time:
                break
            items = copy.deepcopy(my_log_files)
            for log_file in items:
                with open(log_file, 'r') as fd:
                    lines = fd.readlines()
                    if lines and lines[-1].startswith('Task succeeded'):
                        my_log_files.remove(log_file)
                        logger.info('Done: %s', os.path.basename(log_file))
            time.sleep(interval)
        if my_log_files:
            message = 'Logs for dead or unresponsive tasks: {}'\
                      .format([os.path.basename(_) for _ in my_log_files])
            raise RuntimeError(message)


def ir2_hosts():
    remote_hosts = []
    for i in range(1, 11):
        host = f'lsst-dc{i:02}'
        if host in socket.gethostname():
            continue
        remote_hosts.append(host)
    yield np.random.choice(remote_hosts)

def ssh_device_analysis_pool(task_script, device_names, cwd='.', max_time=1800,
                             remote_hosts=None, verbose=False):
    """
    Submit JH tasks on remote nodes.
    """
    cwd = os.path.abspath(cwd)
    setup = os.environ['LCATR_SETUP_SCRIPT']
    task_runner = TaskRunner(task_script, cwd, setup, verbose=verbose)
    num_tasks = len(device_names)

    if remote_hosts is None:
        remote_hosts = ir2_hosts()

    log_files = set()
    with multiprocessing.Pool(processes=num_tasks) as pool:
        outputs = []
        for device_name, remote_host in zip(device_names, remote_hosts):
            log_file = task_runner.log_file(device_name)
            log_files.add(log_file)
            args = remote_host, log_file, device_name
            outputs.append(pool.apply_async(task_runner, args))
        pool.close()
        pool.join()
        [_.get() for _ in outputs]
    task_runner.monitor_tasks(log_files, max_time=max_time)

    return None
