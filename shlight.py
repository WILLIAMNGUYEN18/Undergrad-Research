# importing the required module
from scipy import integrate
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


print("Defining Functions")

#https://scipython.com/book/chapter-8-scipy/examples/visualizing-the-spherical-harmonics/
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html
#http://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html

# scipy.special.sph_harm(m, n, theta, phi) = <ufunc 'sph_harm'>¶
# m = order of harmonic (3?); m <= n
# n = degree of harmonic (???)
# theta between 0, 2pi
# phi between 0, pi
# returns a complex float
#can use .real for real part

#Creating the array meshgrid
#m, l = 2, 3
#phi = np.linspace(0, np.pi, 100)
#theta = np.linspace(0, 2*np.pi, 100)
#phi, theta = np.meshgrid(phi, theta)

#The Cartesian coordinates of the unit sphere
#x = np.sin(phi) * np.cos(theta)
#y = np.sin(phi) * np.sin(theta)
#z = np.cos(phi)

##Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
#fcolors = sph_harm(m, l, theta, phi).real
#fmax, fmin = fcolors.max(), fcolors.min()
#fcolors = (fcolors - fmin)/(fmax - fmin)

# Set the aspect ratio to 1 so our sphere looks spherical
##fig = plt.figure(figsize=plt.figaspect(1.))
#ax = fig.add_subplot(111, projection='3d')
#ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# Turn off the axis planes
#ax.set_axis_off()
#plt.show()

#azimuthal angle is vert xz or yz
#polar is horizontal xy


#rectangle function
def Rect(n):
    if n <= 1:
        return 1.0
    else:
        return 0.0

print("Hardcoding values")


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
        temp = usef(a, b, gfun, hfun, alpha, beta)[0]
        sampleset.append(temp)
        #print('beta: ' + str(beta))
        print('convolution ' + str(temp))
        if (alpha == 0.0): 
            spreadb.append(beta)
        beta += math.pi / bstep
    alpha += math.pi / astep
print('Plotting Graph')


#ax.set_xlabel('Alpha Axis')
#ax.set_ylabel('Beta Axis')
#ax.set_zlabel('Convolution Axis')
#ax.plot_surface(X,Y,Z)
plt.show()