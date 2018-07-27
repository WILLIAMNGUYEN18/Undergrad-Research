msg = "test123"
print(msg)
from scipy import integrate
import math

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



def func(alpha, beta):
    #return math.cos(theta) * Rect((alpha - theta)/w) * Rect((beta - theta)/w)
    #math.cos(theta)
    #need to pull theta somehow
    #print("turn this into necessary function")


