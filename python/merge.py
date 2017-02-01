# Two command-line arguments:  Name of a file containing a list of the single-sensor
# FITS images and the name of the output file

# Merge a set of single CCD images into a (large) FITS file
# assumed to be in S00, SO1, ... order

import astropy.io.fits as fits
import sys

# Read the list of files
io = open(sys.argv[1], "r")
infiles = []
for line in io:
    infiles.append(line.rstrip())
io.close()

hdulist = fits.HDUList()

# Take primary header from the first file in the list
hdu = fits.open(infiles[0])
hdulist.append(hdu[0])

for infile in infiles:
    hdu = fits.open(infile)
    for i in range(1,17):
        hdulist.append(hdu[i])

hdulist.writeto(sys.argv[2])
