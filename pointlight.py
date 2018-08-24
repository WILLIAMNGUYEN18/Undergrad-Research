# cos(theta) - area light model
# 1 graph of cos(theta)
# graphs of cos(theta) - cos(omega)

# importing the required module
from scipy import integrate
from scipy.integrate import quad, dblquad
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.mlab

# example double integral
# f = lambda y, x: x*y**2
# var1 = integrate.dblquad(f, 0, 2, lambda x: 0, lambda x: 1)
# print(var1)
#    (0.6666666666666667, 7.401486830834377e-15)

# example return
# def f(x):
#    return x*x


print("Defining Functions")
#scipy.integrate.dblquad(func, a, b, gfun, hfun)
#func(y, x) from x = a..b and y = gfun(x)..hfun(x).
# temp = usef(a, b, gfun, hfun, alpha, beta)[0]
#a - b = 0 - math.pi /2
#gfun - hfun =  0 - 2 * math.pi

#integral of phi from gfun (0) - hfun(2pi)
#internal integral of theta from a (0) - b (pi/2)


#integrate over cos(theta)d(theta)d(phi)
#phi has much greater effect
#cos(theta) gives more weight

#math.cos(w) = sin(theta) * sin(alpha) * cos(phi) * cos(beta) + sin(theta) * sin(phi) * sin(alpha) * sin(beta) + cos(theta) * cos(alpha)
#function of equation inside double integral
def func (phi, theta, alpha, beta):   
    return math.cos(theta)

#multiplying by scopetheta here will not be effective.
#function to calculate convolution
def usef(a, b, gfun, hfun, alpha, beta):
    return integrate.dblquad(func, a, b,gfun, hfun, args = (alpha, beta))

print("Hardcoding values")
#range for theta
a = 0.0
#b = 2.0 * math.pi
b = math.pi / 2.0
#range for phi
#def gfun ():
#    return 0
gfun = 0.0

#def hfun ():
#    return math.pi
#hfun =  math.pi / 2.0
hfun = 2.0 * math.pi
#arbitrarily set w
w = math.pi / 24.0
#w = 1
Area = (2 * np.pi * (1 - math.cos(w)))

#setting here initially so that they don't need to be passed to functions
step = math.pi/12.0
sampleset = []
spreada = []
spreadb = []
print('Calculating convolutions over alpha and beta')
#will need  to check if alpha/beta looping her works, or if necessary to do while loop
#separate alpha and beta from the iterators and use integer iterators

astep = 24.0
bstep = 24.0
alpha = 0.0

while alpha <= math.pi / 2.0:
    spreada.append(alpha)
    #print('alpha: ' + str(alpha))
    beta = 0.0
    #beta reset for each alpha
    while beta <= (2.0 * math.pi):
        #need to figure out how to work with 4-variable function
        temp = math.cos(alpha)
        sampleset.append(temp)
        #print('beta: ' + str(beta))
        print('convolution ' + str(temp))
        if (alpha == 0.0): 
            spreadb.append(beta)
        beta += math.pi / bstep
    alpha += math.pi / astep
print('Plotting Graph')

print ('premesh x: ' + str(len(spreada)))
print ('premesh y: ' + str(len(spreadb)))
print ('premesh z: ' + str(len(sampleset)))
#X,Y,Z = np.meshgrid(spreada, spreadb, sampleset)
print (len(sampleset))

#Making matrices from spreada and spreadb arrays
#effectively making the grid of x and y (bringing y to 49)
X,Y = np.meshgrid(spreada, spreadb)
print (X.shape)
print (Y.shape)
#reshaping the sampleset to match the grid arrangment of spreada (x) and spreadb(y)
Z = np.reshape(np.array(sampleset), (len(spreada), len(spreadb)))
#transposing z to match the alpha and beta iteration.
Z = np.transpose(Z)
#Z 2d array such that 
print (Z.shape)

print ('mesh x: ' + str(len(X)))
print ('mesh y: ' + str(len(Y)))
print ('mesh z: ' + str(len(sampleset)))
fig = plt.figure()
ax = fig.gca(projection='3d')


#change alpha to azimuth
ax.set_xlabel('Azimuth Axis')

#change beta to zenith
ax.set_ylabel('Zenith Axis')
#lighting intensity due to
ax.set_zlabel('Convolution Axis')

# Customize the z axis.
#ax.set_zlim(-1, 1)
#need to set same graph scale
ax.plot_surface(X,Y,Z)
plt.show()

