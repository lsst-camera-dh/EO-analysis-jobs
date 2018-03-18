"""
Function to copy eotest results files from fe55 analysis jobs for
a specified EO analysis run.
"""
from __future__ import print_function
import os
import shutil
import siteUtils
from DataCatalog import DataCatalog

__all__ = ['copy_fe55_gain_results']

def copy_fe55_gain_results(fe55_run_number, sensor_id, dest='.',
                           job_name='fe55_raft_analysis'):
    """
    Copy eotest results files with Fe55 gain measurements
    for the desired run to the specified destination directory.

    Parameters
    ----------
    fe55_run_number: str
        eTraveler run number of the EO acquisition.
    sensor_id: str
        LSST_NUM of the sensor for the corresponding eotest results file.
    dest: str ["."]
        Destination directory for copied files.
    job_name: str ['fe55_raft_analysis']
        Harnessed job name that produced the eotest results files.

    Raises
    ------
    RuntimeError:  Raised if not exactly one eotest results file is returned.
    """
    site = siteUtils.getSiteName()

    mirror_type = 'prod' if 'Prod' in os.environ['LCATR_LIMS_URL'] else 'test'

    site_folder = "{0}-{1}/{1}".format(site, mirror_type)

    folder = os.path.join("/LSST/mirror", site_folder,
                          siteUtils.getUnitType(), siteUtils.getUnitId(),
                          str(fe55_run_number), job_name)

    datacat = DataCatalog(folder=folder, site=site)
    query = ' && '.join(('TESTTYPE=="FE55"',
                         'LSST_NUM=="%s"' % sensor_id
                         'LsstId=="%s"' % siteUtils.getUnitId()))
    datasets = datacat.find_datasets(query)
    pattern = '%s_eotest_results.fits' % sensor_id

    file_paths = [item for item in datasets.full_paths()
                  if item.endswith(pattern)]
    if len(file_paths) != 1:
        raise RuntimeError("Exactly 1 eotest results file should be "
                           "found for run %s" % fe55_run_number)
    print("copying Fe55 eotest results file to %s :" % dest)
    for item in file_paths:
        print("  ", item)
        shutil.copy(item, '.')
    return [os.path.basename(x) for x in file_paths]
