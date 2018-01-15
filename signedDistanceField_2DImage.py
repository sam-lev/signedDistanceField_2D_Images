#! /usr/bin/python                                                                                                           

import sys
import os
import struct
import math
from scipy.optimize import curve_fit

from PIL import Image
import rawpy
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

from skfmm import distance as sdf
import numpy as np
from scipy.stats import norm

import argparse

parser = argparse.ArgumentParser(description='Takes image and intensity isosurface')
parser.add_argument('-iso', type=float, nargs = 1, help = 'isosurface value for signed distance field', required = False)
parser.add_argument('-im', type=str, nargs = 1, help ='image path to compute signed distance field over', required = True)
args = parser.parse_args()


imPath = args.im[0]


image = Image.open(imPath)

# take only first channel. image can be gray, rgb, rgba
img = image.split()[0]


# Fit a gaussian to the red channel data
# find the mean of the fitted Gaussian
imgAr = np.array(img)
mu, sig = norm.fit(img)
print("mean: "+str(mu))

# Map each red channel value by the fitted
# Gaussian mean, making the mean value
# the 0 iso-surface.
if args.iso:
    mu = args.iso[0]
imgIsoMap = [ i - mu for i in imgAr]

# Compute signed distance field from values mapped for iso as 0
# using the fast marching method.
sdfImage = sdf(imgIsoMap)
imgSignedDist = Image.fromarray(sdfImage)

newImg = Image.new('P', img.size)
newImg.paste(imgSignedDist, (0,0))
newImg.save(imPath)

