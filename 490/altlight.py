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

#need to clamp to?
#max(func, 0)?
def func (phi, theta, alpha, beta):   
    cosgam = (math.sin(theta) * math.sin(alpha) * math.cos(phi) * math.cos(beta) 
    + math.sin(theta) * math.sin(phi) * math.sin(alpha) * math.sin(beta) 
    + math.cos(theta) * math.cos(alpha)) * Rect(theta/w) * math.sin(theta)

    if(cosgam < 0):
        print("NEGATIVE VALUE: " + cosgam)
    return max(cosgam, 0)

#multiplying by scopetheta here will not be effective.
#function to calculate convolution
def usef(a, b, gfun, hfun, alpha, beta):
    return integrate.dblquad(func, a, b,gfun, hfun, args = (alpha, beta))


#rectangle function
def Rect(n):
    if n <= 1:
        return 1.0
    else:
        return 0.0

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
    beta = 0
    #beta reset for each alpha
    while beta <= (2.0 * math.pi):
        #need to figure out how to work with 4-variable function
        temp = usef(a, b, gfun, hfun, alpha, beta)[0] * 1/Area
        sampleset.append(temp)
        #print('beta: ' + str(beta))
        #print('convolution ' + str(temp))
        if (alpha == 0.0): 
            spreadb.append(beta)
        beta += math.pi / bstep
    alpha += math.pi / astep
print('Plotting Graph')

#beta = 0
#while beta <= (2.0 * math.pi):
#    spreadb.append(beta)
#    alpha = 0
#    while(alpha <= math.pi / 2.0):
#        temp = usef(a, b, gfun, hfun, alpha, beta)[0]
#        sampleset.append(temp)
#        #print('beta: ' + str(beta))
#        print('convolution ' + str(temp))
#        if (beta == 0.0):
#            spreada.append(alpha)
#        alpha += math.pi / astep
#    beta += math.pi / bstep 


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

#Looking at the default lengths of these lists
#x is 12 from pi/2 incremented by pi/24 ()
#y is 49 from 2pi incremented by pi/24 (inclusive)
#z is the sample set of across x and y (12 * 49 = 588)
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

#xi = np.linspace(0, math.pi, num = 10)
#yi = np.linspace(0, 2*math.pi, num = 22)
#zi = matplotlib.mlab.griddata(X, Y, Z, xi, yi, interp='linear')
#https://stackoverflow.com/questions/21132758/scipy-interpolation-valueerror-x-and-y-arrays-must-be-equal-in-length-along-int
#not 1D arrays
#np.asarray(a).shape
#squeeze()
#X.flatten()
#Y.flatten()

#print(np.asarray(X).shape)
#print(np.asarray(Y).shape)
#print(np.asarray(np.asarray(sampleset).flatten(order = 'F')))

#zi = matplotlib.mlab.griddata(X, Y, np.asarray(sampleset).flatten, xi, yi, interp='linear')
#zi = matplotlib.mlab.griddata(xi, yi, np.asarray(sampleset).flatten, xi, yi, interp='linear')

#print (str(len(zi)))

#https://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
#need different scale than 4,3
#ax.plot_surface(X.reshape(4,3), Y.reshape(4,3), Z.reshape(4,3))
#ax.plot_surface(X.reshape(4,3), Y.reshape(4,3), np.asarray(sampleset).reshape(4,3))
#ax.plot_surface(X, Y, zi)
#ax.plot_surface(xi, yi, zi)

#change alpha to azimuth
ax.set_xlabel('Azimuth Axis')

#change beta to zenith
ax.set_ylabel('Zenith Axis')
#lighting intensity due to
ax.set_zlabel('Convolution Axis')
ax.set_zlim(0,1)

# Customize the z axis.
#ax.set_zlim(-1, 1)
#need to set same graph scale
ax.plot_surface(X,Y,Z)
plt.show()



#2pi r h
#r 1
#h = (1 - cos)
#2pi * (1 - cos(theta)) = A
#scale each integral by 1/ A