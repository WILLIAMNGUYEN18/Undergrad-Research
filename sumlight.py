from scipy import integrate
from scipy.integrate import quad, dblquad
import math
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

def Rect(n):
    if n < 1/2 and n > -1/2:
        return 1
    else:
        return 0


#might not do alpha beta, may have to use the 2 range variables used in double variable.
def test1(alpha, beta):
    print ("fuck")
    #return math.cos(theta) * Rect((alpha - theta)/w) * Rect((beta - theta)/w)
    #math.cos(theta)
    #need to pull theta somehow
    #print("turn this into necessary function")

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
beta = 0

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

step = math.pi/12
sampleset = list()
#will need  to check if alpha/beta looping her works, or if necessary to do while loop
for alpha in range(0, math.pi, step):
    for beta in range ( 0, 2* math.pi, step):
        #need to figure out how to work with 4-variable function
        temp = usef( a, b, gfun, hfun)
        print(temp)
        sampleset.append(temp)
        





