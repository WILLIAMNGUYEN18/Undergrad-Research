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
w = np.pi / 24.0

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
#potentially forgetting the weight (sqrt(4pi / 2l + 1))
def overall(l,m,w, theta, phi):
    cosstore = normalization(l) * eulers(m, phi) * lpmv(m,l, math.cos(theta))
    #may not need to multiply by math.pi again
    rectstore = normalization(l) * eulers(m, phi) * lpmv(m,l, Rect(theta))
    return 2 * math.pi * math.sqrt((4 * np.pi) / (2 * l + 1)) * cosstore * rectstore



#theta is an array (linspace) atm. Need random theta values
#lpmv produces a single value. This can hopefully be used as a float    
sphcoscalc = []
overl = []
newl = []
step = np.pi / 22

while l < 3:
    overtheta = [0,0,0]
    newset = []

    angphi = 0
    angtheta = 0
    sphcoscalc = []
    while angtheta <= math.pi / 2:
#        print("l = : " + str(l))
#        print("phi = : " + str(angphi))
#        print("theta = : " + str(angtheta))
        store1 = lpmv(m, l, math.cos(angtheta))
        store2 = lpmv(m, l, Rect(angtheta))
#        print("Cosine Legendre function = : " + str(store1.tolist()))
#        print("Rect Legendre function = : " + str(store2.tolist()))

        sphcoscalc.append( normalization(l) * eulers(m, angphi) * lpmv(m,l, math.cos(angtheta)) )            
        #list index out of range
        over1 = overall(l,m,w,angtheta, angphi)
        #print("Cosine and Rectangle Spherical Harmonic Result: " + str(over1))
        overtheta.append(over1)
        #B9 --> coeffcients (flm * kl0 *pl0(cos(theta)))
        presum = over1 * (normalization(l) * lpmv(m,l, math.cos(angtheta)))
        print("New presummed value inserted into set at l:" + str(l) + " and theta: " + str(angtheta) + str(presum))
        newset.append(presum)
        angphi += step
        angtheta += step
    print("")
    for a in overtheta:
        print(a)

    print("")
    overl.append(overtheta)
    newl.append(newset)
    l += 1


print("")
print("")
print("")
print("checking newl values")
for n in newl:
    print( "list at different l: " + str(n))
    for m in n:
        print(m)

fphitheta = []
for counter1, counter2, counter3 in zip(newl[0],newl[1], newl[2]):
    print("value of at l = 0: " + str(counter1))
    print("value of at l = 1: " + str(counter2))
    print("value of at l = 2: " + str(counter3))
    print("value summed across l: " + str(counter1 + counter2 + counter3))
    fphitheta.append((counter1 + counter2 + counter3))


#graph could be 2d, only increasing for theta. No need for 3D



    #for newset B9 calculations
    #potentially ignoring sum m = -l --> l
    #sum of l= 0 --> infinity might go to the chosen band (3?)
    #can just use sum() function along with iterable (list in this case)
    #already separated into lists for each l. able to calculate for entire theta (phi doesn't matter)

        
    #CHANGE STEP TO SEPARATE VARIABLE THAN W
    #add three newl list values

    #print ("length of coefficients is: " + str(len(newl[1])))
    #count = 0
    #while count < len(newl[1]):
    #    newl[1]
    #    count += 1


    #wrong


    #for b in overl:
    #    print("sph coeffcients are: " + str(b)) 
    #comparing spherical harmonic values to my result
    #sphcosharm = sph_harm(m, l,  phi, theta).real
#    print("Our calculated spherical harmonic of cos(theta) is: " + str(sphcoscalc))
#    print("The pre-computed spherical harmonic of cos(theta) is:" + str(sphcosharm))

#need to reverse the loop to separate out the theta values, while grouping the l values
#alternatively, find a way to parse/organize values



    #QUASHED
    #Only getting 12 calculated values for the 24 precomputed. 
    #iterating over angphi would result in n^2 rather than double
    #cut iterations of linspace from 24 to 12.

    #

    #print ("length of calculated: " + str(len(sphcoscalc)))
    #print ("length of precomputed: " + str(len(sphcosharm)))
    #for i in range(len(sphcoscalc)):
    #    if sphcoscalc[i] == sphcosharm[i]:
    #        print("MATCHING CALCULATED & PRECOMPUTED AT: " + str(i))
    #    else:
    #        print("NOT MATCHING CALCULATING & PRECOMPUTED: " + str(i))
#checking steps of sph_harm
#    for x in np.nditer(theta):
#        print(x)
#    for y in np.nditer(phi):
#        print (y)
    #print("At L = " + str(l) + " The product of Rect and Cos Coefficients are: ")
    #for a in overtheta:
    #    print(a)

    #overl.append(overtheta)
    #print("Across All L, All Coefficients: ")
    #for b in overl:
    #    print(b)
    

#for all coefficients, multiply by Kl0Plm(cosϴ)
#for a in overl:
##    for b in a:
#        newset.append[b * normalization(l)]

#will need to think on how to retain connection with phi and theta for graphing purposes

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
