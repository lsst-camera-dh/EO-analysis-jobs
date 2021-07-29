"""
Module of tools to handle BOT data.
"""
from collections import defaultdict
from astropy.io import fits

def most_common_dark_files(dark_files):
    """
    Find dark files with the most common exposure time and return
    the list of filenames.
    """
    def exptime_value(fits_file):
        """Return the value of the EXPTIME keyword in the PHDU."""
        with fits.open(fits_file) as hdus:
            # Use EXPTIME as the key since those seem to be
            # setpoints and not measured values.
            return hdus[0].header['EXPTIME']

    # Collect files by exposure time.
    files_per_exptime = defaultdict(list)
    for item in dark_files:
        files_per_exptime[exptime_value(item)].append(item)

    nmax = 0
    for exptime, files in files_per_exptime.items():
        if len(files) > nmax:
            nmax = len(files)
            candidate = exptime

    print(exptime_value(files_per_exptime[candidate][0]))
    return files_per_exptime[candidate]

if __name__ == '__main__':
    import glob
    dark_files = glob.glob('/gpfs/slac/lsst/fs3/g/data/jobHarness/jh_stage-test/LCA-10134_Cryostat/LCA-10134_Cryostat-0001/7018D/BOT_acq/v0/50493/dark_dark_*/*R22_S11.fits')
    print(len(dark_files))
    darks = most_common_dark_files(dark_files)
    print(len(darks), darks)
