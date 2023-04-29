"""
===============================================================================
Kerr metric


===============================================================================
- Event horizon at r=2M
- ISCO at r = 6M

@author: Eduard Larrañaga - 2023
===============================================================================
"""

from numpy import sin, cos, sqrt

class BlackHole:
    '''
    Definition of the Black Hole described by Schwarzschild metric
    '''
    def __init__(self, M, a):
        self.M = M
        self.a = a
        self.EH = self.M + sqrt(self.M**2 - self.a**2)
        Z1 = self.M + (self.M**2 - self.a**2)**(1/3)*((self.M + self.a)**(1/3) + (self.M - self.a)**(1/3))
        Z2 = sqrt(3*self.a**2 + Z1**2)
        self.ISCOco = 3*self.M + Z2 - sqrt((3*self.M - Z1)*(3*self.M + Z1 + 2*Z2)) 
        self.ISCOcounter = 3*self.M + Z2 + sqrt((3*self.M - Z1)*(3*self.M + Z1 + 2*Z2))

    def metric(self,x):
        '''
        This procedure contains the Schwarzschild metric non-zero components in 
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
        '''
        This function contains the geodesic equations in Hamiltonian form for 
        the Schwarzschild metric
        ===========================================================================
        Coordinates and momentum components
        t = q[0]
        r = q[1]
        theta = q[2]
        phi = q[3]
        k_t = q[4]
        k_r = q[5]
        k_th = q[6]
        k_phi = q[7]
        ===========================================================================
        Conserved Quantities
        E = - k_t
        L = k_phi
        ===========================================================================
        '''

        # Auxiliar Functions
        Sigma = q[1]**2 + self.a**2 * cos(q[2])**2
        Delta = q[1]**2 - 2*self.M*q[1] + self.a**2

        W = -q[4]*(q[1]**2 + self.a**2) - self.a*q[7] 
        partXi = q[1]**2 + (q[7] + self.a*q[4])**2 + self.a**2 *(1 + q[4]**2)*cos(q[2])**2 + (q[7]*cos(q[2])/sin(q[2]))**2
        Xi = W**2 - Delta*partXi

        dXidE = 2*W*(q[1]**2 + self.a**2) + 2.*self.a*Delta*(q[7] + self.a*q[4]*sin(q[2])**2)
        dXidL = -2*self.a*W - 2*self.a*q[4]*Delta - 2*q[7]*Delta/(sin(q[2])**2)

        dXidr = -4*q[1]*q[4]*W - 2*(q[1] - self.M)*partXi - 2*q[1]*Delta 

        dAdr = (q[1] - self.M)/Sigma - (q[1]*Delta)/(Sigma**2)
        dBdr = -q[1]/Sigma**2
        dCdr = dXidr/(2*Delta*Sigma) - (Xi*(q[1]-self.M))/(Sigma*Delta**2) - q[1]*Xi/(Delta*Sigma**2)

        auxth = (self.a**2) * cos(q[2])*sin(q[2])

        dAdth = Delta*auxth/(Sigma**2)
        dBdth = auxth/(Sigma**2)
        dCdth = ((1+q[4]**2)*auxth + q[7]**2 * cos(q[2])/(sin(q[2])**3) )/Sigma + (Xi/(Delta*Sigma**2))*auxth

        # Geodesics differential equations 
        dtdlmbda = dXidE/(2.*Delta*Sigma)
        drdlmbda = (Delta/Sigma)*q[5]
        dthdlmbda = q[6]/Sigma
        dphidlmbda = - dXidL/(2.*Delta*Sigma)
        
        dk_tdlmbda = 0.
        dk_rdlmbda = -dAdr*q[5]**2 - dBdr*q[6]**2 + dCdr 
        dk_thdlmbda = -dAdth*q[5]**2 - dBdth*q[6]**2 + dCdth 
        dk_phidlmbda = 0.
        
        return [dtdlmbda, drdlmbda, dthdlmbda, dphidlmbda, 
                dk_tdlmbda, dk_rdlmbda, dk_thdlmbda, dk_phidlmbda]





###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')