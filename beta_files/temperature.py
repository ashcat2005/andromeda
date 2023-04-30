import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, sqrt, arccos, pi, log

M = 1
a= 0.


Z1 = M + (M**2 - a**2)**(1/3)*((M + a)**(1/3) + (M - a)**(1/3))
Z2 = sqrt(3*a**2 + Z1**2)
ISCOco = 3*M + Z2 - sqrt((3*M - Z1)*(3*M + Z1 + 2*Z2)) 
ISCOcounter = 3*M + Z2 + sqrt((3*M - Z1)*(3*M + Z1 + 2*Z2))

x0 = sqrt(ISCOco/M)
x1 = 2*cos((arccos(a/M) - pi)/3 )
x2 = 2*cos((arccos(a/M) + pi)/3 )
x3 = -2*cos((arccos(a/M))/3 )

def f(r):
    x = sqrt(r/M)
    c = 3/(2*(x**4)*(x**3 - 3*x + 2*a/M) )
    t1 = x - x0 - 3*a*log(x/x0)/(2*M)
    t2 = -((3*(x1-a/M)**2)/(x1*(x1-x2)*(x1-x3)))*log((x-x1)/(x0-x1))
    t3 = -((3*(x2-a/M)**2)/(x2*(x2-x1)*(x2-x3)))*log((x-x2)/(x0-x2))
    t4 = -((3*(x3-a/M)**2)/(x3*(x3-x1)*(x3-x2)))*log((x-x3)/(x0-x3))
    return c*(t1 + t2 + t3 + t4)


rr = np.linspace(7,20,1000)
plt.figure()
plt.plot(rr,f(rr))
plt.show()