"""
Use ssh to run task scripts on remote nodes that are connected to
a common file system and keep track of remote tasks via log files.
"""
import os
import copy
import time
import json
import socket
import logging
import subprocess
import multiprocessing
from collections import defaultdict
import numpy as np
import siteUtils

__all__ = ['ssh_device_analysis_pool']

logging.basicConfig(format='%(asctime)s %(name)s: %(message)s')

class Ir2Hosts:
    """
    Iterator class to provide lsst-dc* host names in a cyclic fashion.
    """
    def __init__(self):
        bad_nodes = os.environ.get('LCATR_BAD_NODES', '').split('_')
        self.hosts = []
        for i in range(1, 11):
            host = f'lsst-dc{i:02}'
            if host in socket.gethostname() or host in bad_nodes:
                continue
            self.hosts.append(host)
        self.num_hosts = len(self.hosts)
        self.index = 0

    def __next__(self):
        value = self.hosts[self.index % self.num_hosts]
        self.index += 1
        return value

    def __iter__(self):
        return self


def zero_func():
    """
    Return 0 to be used the default value for a pickleable
    defaultdict.
    """
    return 0


class TaskRunner:
    """
    Class to manage parallel submission of jh command-line tasks
    on remote machines via ssh.
    """
    def __init__(self, script, working_dir, setup, max_retries=1,
                 remote_hosts=None, verbose=False):
        """
        Parameters
        ----------
        script: str
            Executable script for the desired task. Its first command
            line argument must be the device name on which the
            task operates.
        working_dir: str
            Working directory where the task script should run.
        setup: str
            Path to bash script that will be sourced to set up the
            runtime environment for the task script
        max_retries: int [1]
            Maximum number of retries for failed task executions.
        remote_hosts: iterable [None]
            Iterable that provides the remote hosts for each process.
            If None, then the lsst-dc* hosts will be used.
        verbose: bool [False]
            Flag to output additional diagnostic info about the task
            execution.
        """
        self.params = script, working_dir, setup
        self.max_retries = max_retries
        self.remote_hosts = Ir2Hosts() if remote_hosts is None else remote_hosts
        self.verbose = verbose
        self.log_dir = os.path.join(working_dir, 'logging')
        os.makedirs(self.log_dir, exist_ok=True)
        self.lcatr_envs = siteUtils.get_lcatr_envs()
        self.task_ids = dict()
        self.log_files = dict()
        self.retries = defaultdict(zero_func)
        self.host_map = None

    def make_log_file(self, task_id, clean_up=True, params=None):
        """
        Create a log filename from the task name and task_id and
        clean up any existing log files in the logging directory.
        """
        if params is None:
            params = self.params
        script, _, _ = params
        task_name = os.path.basename(script).split('.')[0]
        log_file = os.path.join(self.log_dir, f'{task_name}_{task_id}.log')
        self.task_ids[log_file] = task_id
        self.log_files[task_id] = log_file
        if clean_up and os.path.isfile(log_file):
            os.remove(log_file)
        return log_file

    def launch_script(self, remote_host, task_id, *args, niceness=10,
                      params=None, wait=False):
        """
        Function to launch the script as a remote process via ssh.
        """
        logger = logging.getLogger('TaskRunner.launch_script')
        logger.setLevel(logging.INFO)

        if params is None:
            params = self.params
        script, working_dir, setup = params
        log_file = self.log_files[task_id]
        command = f'ssh {remote_host} '
        command += f'"cd {working_dir}; source {setup}; '
        for key, value in self.lcatr_envs.items():
            command += f'export {key}={value}; '
        command += f'(echo; nice -n {niceness} ipython '
        command += f'--HistoryManager.enabled=False {script} {task_id} '
        command += ' '.join([str(_) for _ in args])
        command += r' && echo Task succeeded on \`hostname\`'
        command += r' || echo Task failed on \`hostname\`)'
        if wait:
            command += f' &>> {log_file}"'
        else:
            command += f' &>> {log_file}&"'
        if self.verbose:
            logger.info(command)
        logger.info('Launching %s on %s', script, remote_host)
        subprocess.check_call(command, shell=True)

    def monitor_tasks(self, max_time=None, interval=1):
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

        Raises
        ------
        RuntimeError:  This will be raised if max_time is reached.
        """
        logger = logging.getLogger('TaskRunner.monitor_tasks')
        logger.setLevel(logging.INFO)

        # Poll log files for completion of each task:
        log_files = list(self.log_files.values())
        t0 = time.time()
        failures = []
        while log_files:
            if max_time is not None and time.time() - t0 > max_time:
                break
            items = copy.deepcopy(log_files)
            to_retry = []
            for log_file in items:
                task_id = self.task_ids[log_file]
                with open(log_file, 'r') as fd:
                    lines = fd.readlines()
                    if lines and lines[-1].startswith('Task succeeded'):
                        log_files.remove(log_file)
                        logger.info('Done: %s', os.path.basename(log_file))
                    if lines and lines[-1].startswith('Task failed'):
                        if self.retries[task_id] >= self.max_retries:
                            log_files.remove(log_file)
                            logger.info('Failed: %s after %d attempt(s)',
                                        os.path.basename(log_file),
                                        self.max_retries + 1)
                            failures.append(task_id)
                        else:
                            to_retry.append(task_id)
                            self.retries[task_id] += 1
            if to_retry:
                logger.info('Retrying tasks for: ')
                for item in to_retry:
                    logger.info('  %s', item)
                self.submit_jobs(to_retry, retry=True)
            time.sleep(interval)
        messages = []
        if log_files:
            messages.append('\n  Unresponsive tasks: {}'\
                            .format([self.task_ids[_] for _ in log_files]))
        if failures:
            messages.append('  Failed tasks after {} retries: {}'\
                            .format(self.max_retries, failures))
        if messages:
            raise RuntimeError('\n'.join(messages))

    def stage_data(self, device_map_file='device_list_map.json'):
        """
        Function to dispatch data staging script to the remote hosts.
        """
        # Make inverse index of host -> list of devices and
        # save as json for the staging script to use.
        device_map = defaultdict(list)
        for device_name, host in self.host_map.items():
            device_map[host].append(device_name)
        with open(device_map_file, 'w') as fd:
            json.dump(dict(device_map), fd)
        # Loop over hosts and launch the copy script on each host.
        copy_script = os.path.join(os.environ['EOANALYSISJOBSDIR'],
                                   'python', 'stage_bot_data.py')
        # Set params to override self.params in self.make_log_file
        # and self.launch_script
        params = (copy_script, *self.params[1:])
        # Loop over hosts and launch staging script.
        with multiprocessing.Pool(processes=len(device_map)) as pool:
            workers = []
            for host in device_map:
                if host not in self.log_files:
                    self.make_log_file(host, params=params)
                args = (host, host)
                kwds = dict(params=params, wait=True)
                time.sleep(0.5)
                workers.append(pool.apply_async(self.launch_script, args, kwds))
            pool.close()
            pool.join()
            _ = [_.get() for _ in workers]

        # Check for staging failures.
        for log_file in self.log_files.values():
            with open(log_file) as fobj:
                lines = fobj.readlines()
                if lines[-1].startswith('Task failed'):
                    raise RuntimeError(f'Data staging failed: {log_file}')

        # Clear self.log_files of staging script entries.
        self.log_files = dict()

    def submit_jobs(self, device_names, retry=False):
        """
        Submit a task script process for each device.

        Parameters
        ----------
        device_names: list
            List of devices for which the task script will be run.
        """
        num_tasks = len(device_names)
        if not retry:
            self.host_map = dict(zip(device_names, self.remote_hosts))
            if os.environ.get('LCATR_STAGE_DATA', 'True') == 'True':
                self.stage_data()

        # Using multiprocessing allows one to launch the scripts much
        # faster since it can be done asynchronously.
        with multiprocessing.Pool(processes=num_tasks) as pool:
            outputs = []
            for device_name, remote_host in self.host_map.items():
                if device_name not in device_names:
                    # This must be a retry and this device is not
                    # in the list of devices to retry, so skip it.
                    continue
                if device_name not in self.log_files:
                    self.make_log_file(device_name)
                args = remote_host, device_name
                time.sleep(0.5)
                outputs.append(pool.apply_async(self.launch_script, args))
            pool.close()
            pool.join()
            _ = [_.get() for _ in outputs]


def ssh_device_analysis_pool(task_script, device_names, cwd='.', setup=None,
                             max_time=None, remote_hosts=None, verbose=False,
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
        variable will be used, if it is set. Otherwise, `$INST_DIR/setup.sh`
        will be used.
    max_time: float [None]
        Maximum execution time in seconds for the parent task. If this
        time is exceeded, then a RuntimeError will be raised.  If None,
        then no time limit is imposed.
    remote_hosts: iterable [None]
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
        setup = os.environ.get('LCATR_SETUP_SCRIPT',
                               os.path.join(os.environ['INST_DIR'], 'setup.sh'))

    task_runner = TaskRunner(task_script, cwd, setup, max_retries=max_retries,
                             remote_hosts=remote_hosts, verbose=verbose)
    # In order to limit memory usage and to avoid overloading the file
    # server, divide into 2 batches by default if processing for more
    # than 100 devices is requested.
    ndev = len(device_names)
    num_batches = 2 if ndev > 100 else 1

    # Use override value from LCATR_NUM_BATCHES if it is set.
    num_batches = int(os.environ.get('LCATR_NUM_BATCHES', num_batches))
    print("# devices, # batches, # hosts:",
          ndev, num_batches, task_runner.remote_hosts.num_hosts)

    bounds = np.linspace(0, ndev, num_batches + 1, dtype=int)
    print(bounds)
    for imin, imax in zip(bounds[:-1], bounds[1:]):
        task_runner.submit_jobs(device_names[imin:imax])
        task_runner.monitor_tasks(max_time=max_time)
