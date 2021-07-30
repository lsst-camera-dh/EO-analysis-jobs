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
