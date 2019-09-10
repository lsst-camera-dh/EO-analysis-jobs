"""
Use ssh to run task scripts on remote nodes that are connected to
a common file system and keep track of remote tasks via log files.
"""
import os
import copy
import time
import socket
import logging
import subprocess
import multiprocessing
from parsl_execution import get_lcatr_envs

__all__ = ['ssh_device_analysis_pool']

logging.basicConfig(format='%(asctime)s %(name)s: %(message)s')


def ir2_hosts():
    """
    Return a generator that serves up lsst-dc* remote hosts in
    sequence, excluding the current host.
    """
    remote_hosts = []
    for i in range(1, 11):
        host = f'lsst-dc{i:02}'
        if host in socket.gethostname():
            continue
        remote_hosts.append(host)
    num_hosts = len(remote_hosts)
    index = -1
    while True:
        index += 1
        yield remote_hosts[index % num_hosts]


class TaskRunner:
    """
    Class to manage parallel submission of jh command-line tasks
    on remote machines via ssh.
    """
    def __init__(self, script, working_dir, setup, verbose=False):
        """
        Parameters
        ----------
        script: str
            Executable script for the desired task. It should take
            one command line argument, the device name on which the
            task operates.
        working_dir: str
            Working directory where the task script should run.
        setup: str
            Path to bash script that will be sourced to set up the
            runtime environment for the task script
        verbose: bool [False]
            Flag to output additional diagnostic info about the task
            execution.
        """
        self.params = script, working_dir, setup
        self.verbose = verbose
        self.log_dir = os.path.join(working_dir, 'logging')
        os.makedirs(self.log_dir, exist_ok=True)
        self.lcatr_envs = get_lcatr_envs()
        self._task_ids = dict()
        self._log_files = dict()
        self.failures = []

    def make_log_file(self, task_id, clean_up=True):
        """
        Create a log filename from the task_id and clean up any
        existing log files in the working directory.
        """
        script, _, _ = self.params
        task_name = os.path.basename(script).split('.')[0]
        my_log_file = os.path.join(self.log_dir, f'{task_name}_{task_id}.log')
        self._task_ids[my_log_file] = task_id
        self._log_files[task_id] = my_log_file
        if clean_up and os.path.isfile(my_log_file):
            os.remove(my_log_file)
        return my_log_file

    def __call__(self, remote_host, *args):
        """
        Function call-back for launching the remote process via ssh.
        """
        logger = logging.getLogger('TaskRunner.__call__')
        logger.setLevel(logging.INFO)

        script, working_dir, setup = self.params
        task_id = args[0]
        log_file = self._log_files[task_id]
        command = f'ssh {remote_host} '
        command += f'"cd {working_dir}; source {setup}; '
        for key, value in self.lcatr_envs.items():
            command += f'export {key}={value}; '
        command += f'({script} '
        command += ' '.join([str(_) for _ in args])
        command += ' && echo Task succeeded || echo Task failed)'
        command += f' >& {log_file}&"'
        if self.verbose:
            logger.info(command)
        logger.info('Launching %s on %s', script, remote_host)
        subprocess.check_call(command, shell=True)

    def monitor_tasks(self, max_time=None, interval=1, logger=None):
        """
        Function to keep track of remote processes by inspecting the
        log files.

        Parameters
        ----------
        max_time: float [None]
            Maximum time allowed for the parent task to complete.
            Note that processes on the remote nodes are not killed
            if the parent task times out.
        interval: float [1]
            Polling interval in seconds for checking log files for
            task completion.
        logger: logging.Logger [None]
            Logger object. If None, then make one.

        Raises
        ------
        RuntimeError:  This will be raised if max_time is reached.
        """
        if logger is None:
            logger = logging.getLogger('TaskRunner.monitor_tasks')
            logger.setLevel(logging.INFO)

        # Poll log files for completion of each task:
        my_log_files = list(self._log_files.values())
        failures = []
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
                    if lines and lines[-1].startswith('Task failed'):
                        my_log_files.remove(log_file)
                        logger.info('Failed: %s', os.path.basename(log_file))
                        failures.append(self._task_ids[log_file])
            time.sleep(interval)
        if my_log_files:
            message = 'Logs for dead or unresponsive tasks: {}'\
                      .format([os.path.basename(_) for _ in my_log_files])
            raise RuntimeError(message)
        return failures

    def submit_jobs(self, device_names, remote_hosts=None):
        """
        Submit a task script process for each device.

        Parameters
        ----------
        device_names: list
            List of devices for which the task script will be run.
        remote_hosts: generator or list
            Iterable providing the remote hosts on which to run each
            process.
        """
        num_tasks = len(device_names)
        if remote_hosts is None:
            remote_hosts = ir2_hosts()

        with multiprocessing.Pool(processes=num_tasks) as pool:
            outputs = []
            for device_name, remote_host in zip(device_names, remote_hosts):
                self.make_log_file(device_name)
                args = remote_host, device_name
                time.sleep(0.1)
                outputs.append(pool.apply_async(self, args))
            pool.close()
            pool.join()
            _ = [_.get() for _ in outputs]


def ssh_device_analysis_pool(task_script, device_names, cwd='.', setup=None,
                             max_time=1800, remote_hosts=None, verbose=False,
                             max_retries=1):
    """
    Submit JH tasks on remote nodes.

    Parameters
    ----------
    task_script: str
        Path to the script that executes the jh task.
    device_names: list
        List of devices on which to run the jh task.
    cwd: str ['.']
        Path to the working directory. This will be expanded to an abspath.
    setup: str [None]
        Setup script for task runtime environment.  If None, then
        the script pointed to by the `LCATR_SETUP_SCRIPT` environment
        variable will be used.
    max_time: float [1800]
        Maximum execution time in seconds for all subprocess. If this time
        is exceeded, then a RuntimeError will be raised.
    remote_hosts: list-like [None]
        Iterable that provides the remote hosts for each process. If None,
        then the lsst-dc* hosts will be used.
    verbose: bool [False]
        Flag for additional diagnostic output.
    max_retries: int [1]
        Maximum number of retries for failed tasks.

    Raises
    ------
    RuntimeError:  This will be raised if max_time is reached.
    """
    cwd = os.path.abspath(cwd)
    if setup is None:
        setup = os.environ['LCATR_SETUP_SCRIPT']

    task_runner = TaskRunner(task_script, cwd, setup, verbose=verbose)
    retries = 0
    devices_todo = copy.deepcopy(device_names)

    while devices_todo and retries <= max_retries:
        task_runner.submit_jobs(devices_todo, remote_hosts=remote_hosts)
        devices_todo = task_runner.monitor_tasks(max_time=max_time)
        retries += 1
