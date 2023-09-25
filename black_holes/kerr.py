"""
===============================================================================
Kerr metric


===============================================================================
- Event horizon at r = M + sqrt(M^2 - a^2)
- ISCO at r = 3M + Z2 -+ sqrt((3M - Z1)(3M + Z1 + 2Z2))

@author: Eduard Larrañaga - 2023
===============================================================================
"""

from numpy import sin, cos, sqrt, ctypeslib, double, zeros, array
from ctypes import *


class BlackHole:
    '''
    Definition of the Black Hole described by Kerr metric
    '''
    def __init__(self, M, a):
        self.M = M
        self.a = a
        self.EH = self.M + sqrt(self.M**2 - self.a**2)
        Z1 = self.M + (self.M**2 - self.a**2)**(1/3)*((self.M + self.a)**(1/3) + (self.M - self.a)**(1/3))
        Z2 = sqrt(3*self.a**2 + Z1**2)
        self.ISCOco = 3*self.M + Z2 - sqrt((3*self.M - Z1)*(3*self.M + Z1 + 2*Z2)) 
        self.ISCOcounter = 3*self.M + Z2 + sqrt((3*self.M - Z1)*(3*self.M + Z1 + 2*Z2))

        #importation of script in C
        self.dso = CDLL('./black_holes/kbh.so')
        arg_1d   = ctypeslib.ndpointer(dtype=double, ndim=1, flags='C_CONTIGUOUS')
        self.dso.kerrBH.argtypes = [arg_1d, arg_1d, c_double, c_double]


    def metric(self,x):
        '''
        This procedure contains the Kerr metric non-zero components in 
        spherical coordinates
        ===========================================================================
        Coordinates 
        t = x[0]
        r = x[1]
        theta = x[2]
        phi = x[3]
        ===========================================================================
        '''
        # Auxiliary functions
        Delta = x[1]**2 - 2*self.M*x[1] + self.a**2
        Sigma = x[1]**2 + (self.a*cos(x[2]))**2
        
        # Metric components
        g_tt = -(1 - 2*self.M*x[1]/Sigma)
        g_rr = Sigma/Delta
        g_thth = Sigma
        g_phph = (x[1]**2 + self.a**2 + 2*(self.a**2)*self.M*x[1]*sin(x[2])**2/Sigma)* sin(x[2])**2
        g_tph = -2*self.a*self.M*x[1]*sin(x[2])**2/Sigma
        
        return [g_tt, g_rr, g_thth, g_phph, g_tph]

    def geodesics(self, q, lmbda):

        # Auxiliar Functions
        result = zeros(8)
        self.dso.kerrBH(array(q), result, self.M, self.a)
        
        return result.tolist()



###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')