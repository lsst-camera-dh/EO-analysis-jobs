#!/usr/bin/env ipython
"""
Script for BOT nonlinearity analysis.
"""
def nonlinearity_jh_task(det_name):
    """JH version of single sensor execution of the non-linearity task."""
    import glob
    import siteUtils
    from bot_eo_analyses import make_file_prefix, nonlinearity_task

    run = siteUtils.getRunNumber()
    file_prefix = make_file_prefix(run, det_name)

    try:
        detresp_file \
            = siteUtils.dependency_glob(f'{det_name}*_det_response.fits',
                                        jobname='flat_pairs_BOT')[0]
    except IndexError:
        print("nonlinearity_task: detector response file not found for",
              det_name)
        return None

    outfile = f'{file_prefix}_nlc.fits'
    return nonlinearity_task(run, det_name, detresp_file, outfile)

if __name__ == '__main__':
    import sys
    det_name = sys.argv[1]
    nonlinearity_jh_task(det_name)
