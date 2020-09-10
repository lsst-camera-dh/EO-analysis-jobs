"""
Functions to process raft-level single sensor flats for tearing detection.
"""
import os
import matplotlib.pyplot as plt
import lsst.eotest.sensor as sensorTest
import lsst.eotest.image_utils as imutils
import siteUtils

__all__ = ['tearing_detection', 'persist_tearing_png_files']

def tearing_detection(fitsfiles, num_png_files=1, bias_frame=None):
    """
    Run the tearing detection code over a collection of single
    sensor flats.

    Parameters
    ----------
    fitsfiles: list
        List of single sensor flats.  This should generally be for a
        given sensor since any png files of the tearing profiles would
        be later persisted by slot number.
    num_png_files: int [1]
        Number of the tearing profiles to output as a png file.
    bias_frame: str [None]
        Name of bias frame file. If None, then the overscan region will
        be used.

    Returns
    -------
    tuple(list, list, dict): The first item is a list of files for
    which tearing has been detected.  The second item is a list of png
    files of the tearing profiles for the sensors with tearing, if
    requested.  The third item is a dictionary of the number of
    tearing detections per amp, keyed by amp.
    """
    amp_counts = {amp: 0 for amp in imutils.allAmps(fits_file=fitsfiles[0])}

    files_with_tearing = []
    png_files = []
    for filename in fitsfiles:
        ts = sensorTest.TearingStats(filename, bias_frame=bias_frame)
        for amp in amp_counts:
            amp_counts[amp] = max(ts.amp_tearing_count(amp), amp_counts[amp])
        if ts.has_tearing():
            files_with_tearing.append(filename)
            if len(png_files) < num_png_files:
                ts.plot_profiles()
                outfile = (os.path.basename(filename).split('.')[0] +
                           '_tearing.png')
                plt.savefig(outfile)
                png_files.append(outfile)
    return files_with_tearing, png_files, amp_counts


def persist_tearing_png_files(png_files, folder=None, metadata=None):
    """
    Create the lcatr.schema.filerefs for persisting the png files
    with the DataCatalog.

    Parameters
    ----------
    png_files: list
        A list of png files for the tearing profile plots.
    folder: str [None]
        Folder under which to persist the png file.  For raft-level
        analysis, this would be the slot number of the CCD.
    metadata: dict [None]
        Any additional metadata to persist with the png files.

    Returns
    -------
    list: This is a list of lcatr.schema.filerefs.
    """
    if metadata is None:
        metadata = dict()
    md = siteUtils.DataCatalogMetadata(**metadata)
    png_filerefs = []
    for png_file in png_files:
        dp = 'tearing_profiles'
        lsst_id = os.path.basename(png_file).split('_')[0]
        png_filerefs.append(siteUtils.make_fileref(png_file, folder=folder,
                                                   metadata=md(DATA_PRODUCT=dp,
                                                               LsstId=lsst_id)))
    return png_filerefs
