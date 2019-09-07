"""
Function to parallelize subcomponent analyses for an assembly such
as raft or full focal plane.
"""
import os
import logging
from parsl.app.app import python_app, bash_app
import siteUtils
import camera_components
from parsl_ir2_dc_config import load_ir2_dc_config, MAX_PARSL_THREADS


__all__ = ['parsl_sensor_analyses', 'parsl_device_analysis_pool']


@bash_app
def bash_wrapper(script_name, *args, cwd=None, lcatr_envs=None, **kwds):
    script_lines = []
    if cwd is not None:
        script_lines.append(f'cd {cwd}')
    if lcatr_envs is not None:
        script_lines.extend(['export {}={}'.format(*_)
                             for _ in lcatr_envs.items()])
    script_lines.append(' '.join([script_name] + list(args)))
    return '\n'.join(script_lines)


@python_app
def python_wrapper(func, *args, cwd=None, lcatr_envs=None, logger=None,
                   walltime=1800, **kwargs):
    """
    Parsl python_app function wrapper that is serialized and executed
    on worker nodes.
    """
    if cwd is not None:
        # cd to the specified working directory which should be set to
        # the "staging area" for the harnessed job.
        os.chdir(cwd)

    if lcatr_envs is not None:
        # Update the enviroment variables with the LCATR_* values.
        # These can change from one job to the next, so they have to
        # be updated here, just prior to running the harnessed job
        # function.
        os.environ.update(lcatr_envs)

    result = func(*args, **kwargs)
    if logger is None:
        logger = logging.getLogger('python_wrapper')
        logger.setLevel(logging.INFO)

    logger.info('pid {} returning'.format(os.getpid()))
    return result


def get_lcatr_envs():
    """
    Extract the LCATR_* environment variables for updating the runtime
    environment of the harnessed job code to be executed in
    parsl_wrapper.
    """
    lcatr_envs = dict()
    for key, value in os.environ.items():
        if key.startswith('LCATR'):
            lcatr_envs[key] = value
    return lcatr_envs


def parsl_device_analysis_pool(task_func, device_names, processes=None,
                               cwd=None, walltime=3600):
    """
    Use a multiprocessing.Pool to run a device-level analysis task
    over a collection of device names.  The task_func should be
    implemented as pickleable function that takes the desired device
    name as its single argument.

    Parameters
    ----------
    task_func: function (or str)
        A pickleable function that takes the detector name string as
        its argument.  If task_func is a string, then it is interpreted
        as the path of the command-line version of the task.
    device_names: list
        The list of device names to run in the pool.
    processes : int [None]
        The maximum number of processes to have running at once.  If
        None, then set to 1 or one less than the number of cores,
        whichever is larger.
    cwd : str [None]
        The working directory to cd to at the remote node.  Nominally, this
        is a location on the shared file system.
    walltime: float [3600]
        Walltime in seconds for python app execution.  If the python app
        does not return within walltime, a parsl.app.errors.AppTimeout
        exception will be thrown.

    Raises
    ------
    parsl.app.errors.AppTimeout
    """
    load_ir2_dc_config()

    if processes is None:
        # Use the maximum number of cores available, reserving one for
        # the parent process.
        processes = max(1, MAX_PARSL_THREADS - 1)
    processes = int(os.environ.get('LCATR_PARALLEL_PROCESSES', processes))

    print("Running in %i processes" % processes)

    if processes == 1:
        # For cases where only one process will be run at a time, it's
        # faster to run serially instead of using a
        # multiprocessing.Pool since the pickling that occurs can
        # cause significant overhead.
        for device_name in device_names:
            task_func(device_name)
        return None

    # Put the AppFutures in a list so that the task_funcs can run
    # asynchronously on the workers.
    parsl_wrapper = bash_wrapper if isinstance(task_func, str) \
                    else python_wrapper
    outputs = []
    for device_name in device_names:
        print(f"launching parsl job for {task_func} and {device_name}")
        outputs.append(parsl_wrapper(task_func, device_name, cwd=cwd,
                                     lcatr_envs=get_lcatr_envs(),
                                     logger=None, walltime=walltime))

    # Check the resolution of the AppFutures by asking for the result.
    # Calling .result() blocks until the function has exited on the worker node.
    return [_.result() for _ in outputs]


def parsl_sensor_analyses(run_task_func, raft_id=None, processes=None,
                          cwd=None, walltime=3600):
    """
    Run a sensor-level analysis task implemented as a pickleable
    function that takes the desired sensor id as its single argument.

    Parameters
    ----------
    run_task_func : function
        A pickleable function that takes the sensor_id string as
        its argument.
    raft_id : str [None]
        The RTM (or RSA) LSST ID.  If None (default), the LCATR_UNIT_ID
        is used.
    processes : int [None]
        The maximum number of processes to have running at once.
        If None (default), then set to 1 or one less than
        the number of cores, whichever is larger.
    cwd : str [None]
        The working directory to cd to at the remote node.  Nominally, this
        is a location on the shared file system.
    walltime: float [3600]
        Walltime in seconds for python app execution.  If the python app
        does not return within walltime, a parsl.app.errors.AppTimeout
        exception will be thrown.

    Raises
    ------
    parsl.app.errors.AppTimeout
    """
    load_ir2_dc_config()

    if raft_id is None:
        raft_id = siteUtils.getUnitId()

    raft = camera_components.Raft.create_from_etrav(raft_id)

    return parsl_device_analysis_pool(run_task_func, raft.sensor_names,
                                      processes=processes, cwd=cwd,
                                      walltime=walltime)
