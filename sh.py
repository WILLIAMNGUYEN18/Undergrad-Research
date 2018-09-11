# importing the required module
from scipy import integrate
from scipy import special
from scipy.special import sph_harm
from scipy.integrate import quad, dblquad
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.mlab
import pyshtools


theta = np.linspace(0, np.pi / 2, 48)
phi = np.linspace(0, 2 * np.pi, 48)

#theta1, phi1 = np.meshgrid(theta, phi)
m = 0
l = 0
w = np.pi / 24.0

sampleset = []

#unsure if I need these
spreada = []
spreadb = []

def normalization(l):
    #2 possibilities here based on readings
    #return
    #norm = math.sqrt((2 * l + 1) / (4 * math.pi))
    norm = math.sqrt((4 * math.pi) / (2 * l + 1))

    return norm

def eulers(m,phi):
    eul = 1
    return eul
    #m always equals 0, so euler's function always equals e^0 == 1
    #ignoring imaginary piece (i) in sample
    #eul = math.pow(math.e, math.com m * phi)

def Rect(n):
    if n <= 1:
        return 1.0
    else:
        return 0.0

def overall(l,m,w, theta, phi,):
    cosstore = 2 * math.pi *  normalization(l) * eulers(m, phi)
    #may not need to multiply by math.pi again
    rectstore = 2 * math.pi * normalization(l) * eulers(m, phi)



#scipy.special.lpmv(m, v, x)
#m : array_like
#    Order (int or float). If passed a float not equal to an integer the function returns NaN.
#v : array_like
#    Degree (float).
#x : array_like
#    Argument (float). Must have |x| <= 1.
#returns:
#pmv : ndarray
#    Value of the associated Legendre function.

#may need to create a full algorithm to intake theta and phi while not using it for calculations until the end
    #float is required for lpmv
    #need to input cos(theta) value to get that value

    #rect function intakes theta value (between 0 and pi/2) and returns whether it will be 1 or zero
        #the value may be too high( |x| <= 1)

#unsure how to convert spherical harmonic function back to time-based function
