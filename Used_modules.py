import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import astropy.units as u
from astropy.io import fits
# from modules import rms
import FITS_tools
from astropy import modeling
import matplotlib.pyplot as plt
from numpy import exp, loadtxt, pi, sqrt
from scipy.optimize import curve_fit
from lmfit import Model
from spectral_cube import SpectralCube
from astroquery.esasky import ESASky
from astroquery.utils import TableList
from astropy.wcs import WCS
from reproject import reproject_interp
from numpy import *
from reproject import reproject_interp
import FITS_tools
from scipy.linalg import fractional_matrix_power
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from fits_align.ident import make_transforms
from fits_align.align import affineremap
from glob import glob
from numpy import shape
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from collections import OrderedDict
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.data import get_pkg_data_filename
import aplpy
from scipy.signal import find_peaks
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
import os
#import cv2
import numpy as np
from scipy.stats import norm
import math
import numpy as np
import scipy as sc
from math import log
from decimal import Decimal
import operator
import matplotlib as mpl
import csv
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import matplotlib.pyplot as pyplot
from pylab import meshgrid, cm, imshow, contour, clabel, colorbar, axis, title, show
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.io.idl import readsav
import matplotlib.colors as colors
import scipy.io as spio
# import idlsave
import pylab as pl
from matplotlib import cm
import numpy as np
from scipy import interpolate
from scipy import ndimage
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
import aplpy
from matplotlib.colors import LogNorm
from scipy.linalg import fractional_matrix_power
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
from astropy import units as u
import seaborn as sns
from scipy.stats import iqr
