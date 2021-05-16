# Uncertainty.py 
# Fiona Luo
# This program finds the uncertainty and rms of a single observation using the 
# corr.fits file from nova.astrometry.net. This uncertainty is used is the LuoOD.py code
# for the Monte Carlo simulation. 

from astropy.io import fits
import numpy as np
import math
from math import pi
import odlib as od

# Return rms given two arrays
# The rms is the root of the mean of the square of the differences between 
# expected (field) and measured (index) values.
def rms(fields, indeces):
    sum = 0
    for i in range(len(fields)):
        diff = (fields[i] - indeces[i])**2
        sum += diff
    sum /= len(fields)
    sum = sum ** 0.5
    return sum

# Input filename for fitsfile. Ex: 
def rmsradec(fitsfile, radians = False):
    table = fits.open(fitsfile)[1].data

    field_ra = table.field_ra
    field_dec = table.field_dec
    index_ra = table.index_ra
    index_dec = table.index_dec
    if (radians):
        return rms(field_ra, index_ra)*pi/180, rms(field_dec, index_dec)*pi/180
    return rms(field_ra, index_ra), rms(field_dec, index_dec)