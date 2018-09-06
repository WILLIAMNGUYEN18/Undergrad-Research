# importing the required module
from scipy import integrate
from scipy import special
from scipy.special import sph_harm

#from scipy.special import sph_harm
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

#associated legendre function?
#transform?
#cos(theta)
# Alm = Clm double integral f(theta) Plm(cos(theta))e^(im(phi)) sin(theta)d(theta)d(phi) 
#== Clm integral(0-->2pi) e^(im(phi)) integral f(theta) Plm(cos(theta)) sin(theta)d(theta)
# @ m = 0, e^(im(phi)) = 1 from 0 --> 2pi = 2pi
# @ non-zero cancels itself out and becomes 0
#Light coefficient  Llm= 2piC Cl0 integral (0 --> omega) Pl0(cos(theta))sin(theta)d(theta), @ m = 0
#0 @ m =/= 0
#Shading COefficient Slm = 2pi Cl0 integral(0 --> pi/2)cos(theta) Pl0(cos(theta))sin(theta)d(theta) @ m = 0
#0 @ m =/= 0

#truncate to l (order) --> 3 rather than infinity

#can use sph_harm to skip steps plugging in specific n, m, phi, and theta
#utilize 2pi and potentially sin(theta)
#sum of l = 0 --> 3 (or infinity) 

#alternatively, try to work with legendre functions and a constant()

print("Defining Functions")

#https://scipython.com/book/chapter-8-scipy/examples/visualizing-the-spherical-harmonics/
#https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html
#http://docs.enthought.com/mayavi/mayavi/auto/example_spherical_harmonics.html
#https://www.youtube.com/watch?v=wuCArO5Iy1A


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

#Will need to adjust theta for my model?
theta = np.linspace(0, np.pi / 2, 48)
#phi may be fine as it becomes a non-factor?
phi = np.linspace(0, 2 * np.pi, 48)
theta1, phi1 = np.meshgrid(theta, phi)

#3rd order --> m = 3; l <= m
m = 0
l = 0
w = np.pi / 24.0
values = []


#for m in range(4):
for l in range(4):
    temp = 2 * np.pi * sph_harm(m, l, phi, theta).real
    print(str(temp))
    values.append(temp)

#        values.append(2 * np.pi * sp.sph_harm(l,m,theta, phi).real)
#        2 * np.pi * scipy.special.sph_harm(l,m,theta, phi).real


#sph_harm(m,l, phi, theta).real


#need to mutliply spherical harmonic of cos(theta) by the spherical harmonic of Rect(theta/w)
#will probably need multiple phi and thetas

#https://keisan.casio.com/has10/SpecExec.cgi?path=08000000.Special%2520Function%252F07001000.Orthogonal%2520polynomial%252F10010100.Spherical%2520harmonics%252Fdefault.xml&charset=utf-8

#phi will be based on how many degrees the vector (x1,y1,z1) deviates from x (or y) on the xy plane
#theta will be based on how many degrees the vector (x1,y1,z1) deviates from the the z vector (0,0,1)

#need to have 2 sets of x,y,z values for the two different vectors.
#this will be done through calculating the x,y,z as theta and phi
#this might be done automatically through iterating through alpha beta

#Do we need to skew corresponding to functions? If not, how are the functions included in the spherical harmonics? If so, how would we skew it?

#take samples over the function?



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
        #temp = usef(a, b, gfun, hfun, alpha, beta)[0]
        #sampleset.append(temp)
        #print('beta: ' + str(beta))
        #print('convolution ' + str(temp))
        if (alpha == 0.0): 
            for l in range(4):
                store = 2 * np.pi * sph_harm(m, l, beta, alpha).real
                print(str(store))
                sampleset.append(store)
            spreadb.append(beta)
        beta += math.pi / bstep
    alpha += math.pi / astep
print('Plotting Graph')
X,Y = np.meshgrid(spreada, spreadb)
Z = np.reshape(np.array(sampleset), (len(spreada), len(spreadb)))

#ax.set_xlabel('Alpha Axis')
#ax.set_ylabel('Beta Axis')
#ax.set_zlabel('Convolution Axis')
#ax.plot_surface(X,Y,Z)
plt.show()

#spherical harmonics for rect and cos functions
#evaluating for
#following thesis chapter
#Don't care about m (whatever controls phi), it's always zero
#Summing for different Ls and
#coefficients for ach of rect and cos
#constant correction based on L
#pairwise multiply
#reconvert back into final answer
#function of theta at that poitn
#sample for different thetas

#write out steps mathematically and get it reviewed.
#remember the 2pi
#make sure to apply constants