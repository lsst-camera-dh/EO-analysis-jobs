"""
Function to parallelize subcomponent analyses for an assembly such
as raft or full focal plane.
"""
import os
import multiprocessing
import siteUtils
import camera_components

__all__ = ['sensor_analyses', 'run_sensor_analysis_pool']

def run_sensor_analysis_pool(task_func, det_names, processes=None):
    """
    Use a multiprocessing.Pool to run a sensor-level analysis task
    over a collection of detector names.  The task_funcs should be
    implemented as pickleable function that takes the desired detector
    name as its single argument.

    Parameters
    ----------
    task_func: function
        A pickleable function that takes the detector name string as
        its argument.
    det_names: list
        The list of detector names to run in the pool.
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
        processes = max(1, multiprocessing.cpu_count() - 1)
    processes = int(os.environ.get('LCATR_PARALLEL_PROCESSES', processes))

    if processes == 1:
        # For cases where only one process will be run at a time, it's
        # faster to run serially instead of using a
        # multiprocessing.Pool since the pickling that occurs can
        # cause significant overhead.
        for det_name in det_names:
            task_func(det_name)
    else:
        with multiprocessing.Pool(processes=processes) as pool:
            results = [pool.apply_async(task_func, (det_name,))
                       for det_name in det_names]
            pool.close()
            pool.join()
            for res in results:
                res.get()

def sensor_analyses(run_task_func, raft_id=None, processes=None):
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
    if raft_id is None:
        raft_id = siteUtils.getUnitId()

    raft = camera_components.Raft.create_from_etrav(raft_id)

    run_sensor_analysis_pool(run_task_func, raft.sensor_names,
                             processes=processes)
