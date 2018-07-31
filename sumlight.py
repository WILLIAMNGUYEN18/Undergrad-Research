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
#msg = "test123"
#print(msg)


# example double integral
#f = lambda y, x: x*y**2
#var1 = integrate.dblquad(f, 0, 2, lambda x: 0, lambda x: 1)
#print(var1)
#    (0.6666666666666667, 7.401486830834377e-15)

# example return
# def f(x):
#    return x*x
print("Defining Functions")
#scipy.integrate.dblquad(func, a, b, gfun, hfun)
#func(y, x) from x = a..b and y = gfun(x)..hfun(x).
#1st parameter func 
#, alpha, beta
def func (phi, theta):
    return math.cos(theta) * Rect(alpha - theta / w) * Rect(beta - theta /w)

#phi, theta, alpha, beta, params removed
def usef(a, b, gfun, hfun):
    return integrate.dblquad(func, a, b,gfun, hfun)

#do we need io?

#plug arbitrary input values for alpha beta currently.

#need to sample (Step sampling for alpha and beta for the above function)
#alpha sampled from 0 - pi, and beta from 0 to 2pi
#after sampling, plot the function (matplotlib?)


def Rect(n):
    if n < 1/2 and n > -1/2:
        return 1
    else:
        return 0


#might not do alpha beta, may have to use the 2 range variables used in double variable.
#def test1(alpha, beta):
#   print ("fuck")
    #return math.cos(theta) * Rect((alpha - theta)/w) * Rect((beta - theta)/w)
    #math.cos(theta)
    #need to pull theta somehow
    #print("turn this into necessary function")

print("Hardcoding values")
#range for theta
a = 0
b = math.pi

#range for phi
#def gfun ():
#    return 0
gfun = 0

#def hfun ():
#    return math.pi
hfun = 2 * math.pi

#arbitrarily set w
w = 1

#setting here initially so that they don't need to be passed to functions
alpha = 0
step = math.pi//12
sampleset = list()
counter = 0
spreada = list()
spreadb = list()
print('Calculating convolutions over alpha and beta')
#will need  to check if alpha/beta looping her works, or if necessary to do while loop
#separate alpha and beta from the iterators and use integer iterators
while alpha < math.pi:
    alpha += math.pi / 12
    spreada.append(alpha)
    print('alpha: ' + str(alpha))

    #beta reset for each alpha
    beta = 0
    while beta < (2 * math.pi):
        #need to figure out how to work with 4-variable function
        temp = usef(a, b, gfun, hfun)
        sampleset.append(temp)
        print('beta: ' + str(beta))
        print('convolution ' + str(temp))
        beta += math.pi / 12
        spreadb.append(beta)

print('Plotting Graph')

#ax.plot3D(xline, yline, zline, 'gray')
#plt.plot(spread, sampleset)
#plt.xlabel('alpha&beta value')
#plt.ylabel('convolution value')
#plt.title('Convolution value across alpha&beta values')

#3d graph default needs vertices.
#https://stackoverflow.com/questions/44355270/plot-3d-in-python-using-three-lists

#ax = plt.axes(projection='3d')
#ax.scatter3D(spreada, spreadb, sampleset,c = sampleset,cmap = 'gray')
#X,Y,Z = np.meshgrid(X,Y,Z)
#https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.meshgrid.html
print ('premesh x: ' + str(len(spreada)))
print ('premesh y: ' + str(len(spreadb)))
print ('premesh z: ' + str(len(sampleset)))
#X,Y,Z = np.meshgrid(spreada, spreadb, sampleset)
X,Y = np.meshgrid(spreada, spreadb)

print ('mesh x: ' + str(len(X)))
print ('mesh y: ' + str(len(Y)))
print ('mesh z: ' + str(len(sampleset)))
fig = plt.figure()
ax = fig.gca(projection='3d')

xi = np.linspace(0, math.pi, num = 10)
yi = np.linspace(0, 2*math.pi, num = 22)
#zi = matplotlib.mlab.griddata(X, Y, Z, xi, yi, interp='linear')
#https://stackoverflow.com/questions/21132758/scipy-interpolation-valueerror-x-and-y-arrays-must-be-equal-in-length-along-int
#not 1D arrays
#np.asarray(a).shape
#squeeze()
X.flatten()
Y.flatten()

print(np.asarray(X).shape)
print(np.asarray(Y).shape)
print(np.asarray(np.asarray(sampleset).flatten()))

zi = matplotlib.mlab.griddata(X, Y, np.asarray(sampleset), xi, yi, interp='linear')
#print (str(len(zi)))

#https://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
#need different scale than 4,3
#ax.plot_surface(X.reshape(4,3), Y.reshape(4,3), Z.reshape(4,3))
#ax.plot_surface(X.reshape(4,3), Y.reshape(4,3), np.asarray(sampleset).reshape(4,3))
#ax.plot_surface(X, Y, zi)
ax.plot_surface(xi, yi, zi)

plt.show()