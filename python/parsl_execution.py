"""
Function to parallelize subcomponent analyses for an assembly such
as raft or full focal plane.
"""
import os
#import siteUtils
import camera_components

import parsl

from parsl.app.app import python_app

from parsl_ir2_dc_config import load_ir2_dc_config,\
    MAX_PARSL_THREADS, PARSL_LOADED


__all__ = ['parsl_sensor_analyses', 'parsl_device_analysis_pool']

@python_app
def parsl_wrapper(func, *args, **kwargs):
    return func(*args, **kwargs)


def parsl_device_analysis_pool(task_func, device_names, processes=None):
    """
    Use a multiprocessing.Pool to run a device-level analysis task
    over a collection of device names.  The task_func should be
    implemented as pickleable function that takes the desired device
    name as its single argument.

    Parameters
    ----------
    task_func: function
        A pickleable function that takes the detector name string as
        its argument.
    device_names: list
        The list of device names to run in the pool.
    processes : int, optional
        The maximum number of processes to have running at once.
        If None (default), then set to 1 or one less than
        the number of cores, whichever is larger.

    Notes
    -----
    Exceptions from subprocesses will be buffered until all of the
    subprocesses have finished running.  If any exceptions are thrown
    and are uncaught, a non-zero exit code will be generated for the
    overall process. Output to stdout or stderr from the subprocesses
    will be interleaved.

    Users can override the default or keyword argument values by setting
    the LCATR_PARALLEL_PROCESSES environment variable.
    """
    if processes is None:
        # Use the maximum number of cores available, reserving one for
        # the parent process.
        processes = max(1, MAX_PARSL_THREADS - 1)
    processes = int(os.environ.get('LCATR_PARALLEL_PROCESSES', processes))

    print ("Running in %i processes" % processes)

    if processes == 1:
        # For cases where only one process will be run at a time, it's
        # faster to run serially instead of using a
        # multiprocessing.Pool since the pickling that occurs can
        # cause significant overhead.
        for device_name in device_names:
            task_func(device_name)
    else:
        print (device_names)
        outputs = [ parsl_wrapper(task_func, device_name).result() for device_name in device_names ]
        print (outputs)
    return outputs


def parsl_sensor_analyses(run_task_func, raft_id=None, processes=None):
    """
    Run a sensor-level analysis task implemented as a pickleable
    function that takes the desired sensor id as its single argument.

    Parameters
    ----------
    run_task_func : function
        A pickleable function that takes the sensor_id string as
        its argument.
    raft_id : str, optional
        The RTM (or RSA) LSST ID.  If None (default), the LCATR_UNIT_ID
        is used.
    processes : int, optional
        The maximum number of processes to have running at once.
        If None (default), then set to 1 or one less than
        the number of cores, whichever is larger.
    """
    if not PARSL_LOADED:
        load_ir2_dc_config()

    if raft_id is None:
        raft_id = siteUtils.getUnitId()

    raft = camera_components.Raft.create_from_etrav(raft_id)

    print(raft.sensor_names)
    parsl_device_analysis_pool(run_task_func, raft.sensor_names,
                               processes=processes)


