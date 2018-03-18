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

def copy_fe55_gain_results(fe55_run_number, dest='.', num_ccds=9,
                           job_name='fe55_raft_analysis'):
    """
    Copy eotest results files with Fe55 gain measurements
    for the desired run to the specified destination directory.

    Parameters
    ----------
    fe55_run_number: str
        eTraveler run number of the EO acquisition.
    dest: str ["."]
        Destination directory for copied files.
    num_ccds: int [9]
        Number of eotest results files to retrieve with a default of 9
        per raft.
    job_name: str ['fe55_raft_analysis']
        Harnessed job name that produced the eotest results files.

    Raises
    ------
    RuntimeError:  if the number of eotest results files does not
        match num_ccds.
    """
    site = siteUtils.getSiteName()

    mirror_type = 'prod' if 'Prod' in os.environ['LCATR_LIMS_URL'] else 'test'

    site_folder = "{0}-{1}/{1}".format(site, mirror_type)

    folder = os.path.join("/LSST/mirror", site_folder,
                          siteUtils.getUnitType(), siteUtils.getUnitId(),
                          str(fe55_run_number), job_name)

    datacat = DataCatalog(folder=folder, site=site)
    query = ' && '.join(('TESTTYPE=="FE55"',
                         'LsstId=="%s"' % siteUtils.getUnitId()))
    datasets = datacat.find_datasets(query)
    pattern = 'eotest_results.fits'

    file_paths = [item for item in datasets.full_paths()
                  if item.endswith(pattern)]
    if len(file_paths) != num_ccds:
        raise RuntimeError("Expected number of Fe55 gain results files "
                           "were not found for run %s" % fe55_run_number)
    print("copying Fe55 eotest results files to %s :" % dest)
    for item in file_paths:
        print("  ", item)
        shutil.copy(item, '.')
    return [os.path.basename(x) for x in file_paths]

if __name__ == '__main__':
    os.environ['SITENAME'] = 'SLAC'
    os.environ['LCATR_LIMS_URL'] = 'Prod'
    os.environ['LCATR_UNIT_ID'] = 'LCA-11021_RTM-002_ETU1'
    os.environ['LCATR_UNIT_TYPE'] = 'LCA-11021_RTM'

    fe55_run_number = 7167
    copy_fe55_gain_results(fe55_run_number)
