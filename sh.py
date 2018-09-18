# importing the required module
from scipy import integrate
from scipy import special
from scipy.special import sph_harm
from scipy.special import lpmv
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


theta = np.linspace(0, np.pi / 2, 12)
phi = np.linspace(0, 2 * np.pi, 12)

#theta1, phi1 = np.meshgrid(theta, phi)
m = 0
l = 0
w = np.pi / 22.0

sampleset = []

#unsure if I need these
spreada = []
spreadb = []

def normalization(l):
    #2 possibilities here based on readings
    #return
    norm = math.sqrt((2 * l + 1) / (4 * math.pi))
    #norm = math.sqrt((4 * math.pi) / (2 * l + 1))

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

#2pi * normal (constant?) * e^im(phi) * Plm(legendre function)
def overall(l,m,w, theta, phi):
    cosstore = 2 * math.pi *  normalization(l) * eulers(m, phi) * lpmv(m,l, math.cos(theta))
    #may not need to multiply by math.pi again
    rectstore = 2 * math.pi * normalization(l) * eulers(m, phi) * lpmv(m,l, Rect(theta))
    return cosstore * rectstore



#theta is an array (linspace) atm. Need random theta values
#lpmv produces a single value. This can hopefully be used as a float    
sphcoscalc = []


while l < 3:
    overtheta = [0,0,0]
    angphi = 0
    angtheta = 0
    sphcoscalc = []
    while angtheta <= math.pi / 2:
        print("l = : " + str(l))
        print("phi = : " + str(angphi))
        print("theta = : " + str(angtheta))
        store1 = lpmv(m, l, math.cos(angtheta))
        store2 = lpmv(m, l, Rect(angtheta))
        print("Cosine Legendre function = : " + str(store1.tolist()))
        print("Rect Legendre function = : " + str(store2.tolist()))
        over1 = overall(l,m,w,angtheta, angphi)
        print("Cosine and Rectangle Spherical Harmonic Result: " + str(over1))
        sphcoscalc.append( normalization(l) * eulers(m, angphi) * lpmv(m,l, math.cos(angtheta)) )            
        #list index out of range
        overtheta[l] += over1
        angphi += w
        angtheta += w

    #comparing spherical harmonic values to my result
    sphcosharm = sph_harm(m, l,  phi, theta).real
    print("Our calculated spherical harmonic of cos(theta) is: " + str(sphcoscalc))
    print("The pre-computed spherical harmonic of cos(theta) is:" + str(sphcosharm))

    #QUASHED
    #Only getting 12 calculated values for the 24 precomputed. 
    #iterating over angphi would result in n^2 rather than double
    #cut iterations of linspace from 24 to 12.

    #

    print ("length of calculated: " + str(len(sphcoscalc)))
    print ("length of precomputed: " + str(len(sphcosharm)))
    for i in range(len(sphcoscalc)):
        if sphcoscalc[i] == sphcosharm[i]:
            print("MATCHING CALCULATED & PRECOMPUTED AT: " + str(i))
        else:
            print("NOT MATCHING CALCULATING & PRECOMPUTED: " + str(i))
    for x in np.nditer(theta):
        print(x)
    
    for y in np.nditer(phi):
        print (y)
    l += 1



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
