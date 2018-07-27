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

#range for psi
gfun = 0
hfun = 2 * math.pi


#scipy.integrate.dblquad(func, a, b, gfun, hfun
#1st parameter func 
def func (theta, psi):
    return math.cos(theta) * Rect(alpha - theta / w) * Rect(beta - theta /w)

integrate.dblquad(func, a, b,gfun, hfun)





