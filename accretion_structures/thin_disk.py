"""
===============================================================================
Novikov-Thorne thin accretion disk with a time-averaged energy flux
emitted from the surface of the disk
===============================================================================
@author: Eduard Larrañaga - 2023
===============================================================================
"""
from numpy import sin, cos, sqrt, arccos, pi, log

class structure:
    def __init__(self, blackhole, R_min , R_max, corotating=True):
        self.in_edge = R_min
        self.out_edge = R_max
        self.M = blackhole.M
        self.a = blackhole.a
        if corotating:
            ISCO = blackhole.ISCOco
        else:
            ISCO = blackhole.ISCOcounter
        self.x0 = sqrt(ISCO/self.M)
        self.x1 = 2*cos((arccos(self.a/self.M) - pi)/3 )
        self.x2 = 2*cos((arccos(self.a/self.M) + pi)/3 )
        self.x3 = -2*cos((arccos(self.a/self.M))/3 )

    def energy_flux(self, r):
        if r>self.in_edge and r<self.out_edge:
            x = sqrt(r/self.M)
            c = 3/(2*(x**4)*(x**3 - 3*x + 2*self.a/self.M) )
            t1 = x - self.x0 - 3*self.a*log(x/self.x0)/(2*self.M)
            t2 = -((3*(self.x1-self.a/self.M)**2)/(self.x1*(self.x1-self.x2)*(self.x1-self.x3)))*log((x-self.x1)/(self.x0-self.x1))
            t3 = -((3*(self.x2-self.a/self.M)**2)/(self.x2*(self.x2-self.x1)*(self.x2-self.x3)))*log((x-self.x2)/(self.x0-self.x2))
            t4 = -((3*(self.x3-self.a/self.M)**2)/(self.x3*(self.x3-self.x1)*(self.x3-self.x2)))*log((x-self.x3)/(self.x0-self.x3))
            return c*(t1 + t2 + t3 + t4)
        else:
            return 0.



###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')
