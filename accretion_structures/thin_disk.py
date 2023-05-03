"""
===============================================================================
Novikov-Thorne thin accretion disk with a time-averaged energy flux
emitted from the surface of the disk
===============================================================================
@author: Eduard Larrañaga - 2023
===============================================================================
"""
from numpy import sin, cos, sqrt, arccos, pi, log, linspace, min
from scipy.interpolate import interp1d

class structure:
    def __init__(self, blackhole, corotating=True, R_min=False, R_max=20.):
        self.out_edge = R_max
        self.M = blackhole.M
        self.a = blackhole.a
        if corotating:
            self.ISCO = blackhole.ISCOco
            self.in_edge = blackhole.ISCOco
        else:
            self.ISCO = blackhole.ISCOcounter
            self.in_edge = blackhole.ISCOcounter
        if R_min:
            self.in_edge = R_min

        rr = linspace(self.in_edge, self.out_edge, 100000)
        ff = self.f(rr)
        ff = ff - min(ff)
        self.energy = interp1d(rr,ff)

    def f(self, r):
        x0 = sqrt(self.ISCO/self.M)
        x1 = 2*cos((arccos(self.a/self.M) - pi)/3 )
        x2 = 2*cos((arccos(self.a/self.M) + pi)/3 )
        x3 = -2*cos((arccos(self.a/self.M))/3 )
        x = sqrt(r/self.M)
        c = 3/(2*(x**4)*(x**3 - 3*x + 2*self.a/self.M) )
        t1 = x - x0 - 3*self.a*log(x/x0)/(2*self.M)
        t2 = -((3*(x1-self.a/self.M)**2)/(x1*(x1-x2)*(x1-x3)))*log((x-x1)/(x0-x1))
        t3 = -((3*(x2-self.a/self.M)**2)/(x2*(x2-x1)*(x2-x3)))*log((x-x2)/(x0-x2))
        t4 = -((3*(x3-self.a/self.M)**2)/(x3*(x3-x1)*(x3-x2)))*log((x-x3)/(x0-x3))
        return c*(t1 + t2 + t3 + t4)

    
    def energy_flux(self,r):
        if r>self.in_edge and r<self.out_edge:
            return self.energy(r)
        else:
            return 0.




###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')
