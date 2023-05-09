"""
===============================================================================
Scalar Hair Black Hole 

ds^2 = g_tt dt^2 + g_rr dr^2  + r^2 dtheta^2 + r^2 sin^2 (theta) dphi^2

===============================================================================
- Event horizon at r = 1

@author: Eduard Larrañaga - 2023
===============================================================================
"""

from numpy import sin, cos, loadtxt, linspace, asarray
from scipy.interpolate import interp1d

class BlackHole:
    '''
    Definition of the Black Hole described by Schwarzschild metric
    '''
    def __init__(self, M):
        data = loadtxt('metrics/numerical_data/scalarBH/metricpp0=0.1.txt')
        self.g_tt = interp1d(data[:,0], data[:,1], bounds_error=False, fill_value = 0)
        self.g_rr = interp1d(data[:,0], data[:,2], bounds_error=False, fill_value = 0)
        self.gtt = interp1d(data[:,0], data[:,3], bounds_error=False, fill_value = 0)
        self.grr = interp1d(data[:,0], data[:,4], bounds_error=False, fill_value = 0)
        self.drgtt = interp1d(data[:,0], data[:,5], bounds_error=False, fill_value = 0)
        self.drgrr = interp1d(data[:,0], data[:,6], bounds_error=False, fill_value = 0)
        self.M = M
        self.a = 0.
        self.EH = 2*M
        self.ISCOco = 6*M 
        self.ISCOcounter = 6*M

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
        # Metric components
        g_tt = self.g_tt(x[1])
        g_rr = self.g_rr(x[1])
        g_thth = x[1]**2
        g_phph = (x[1]*sin(x[2]))**2
        g_tph = 0.
        
        return [g_tt, g_rr, g_thth, g_phph, g_tph]
    
    def inverse_metric(self,x):
        '''
        This procedure contains the Schwarzschild inverse metric non-zero 
        components in spherical coordinates
        ===========================================================================
        Coordinates 
        t = x[0]
        r = x[1]
        theta = x[2]
        phi = x[3]
        ===========================================================================
        '''
        # Metric components
        gtt = self.gtt(x[1])
        grr = self.grr(x[1])
        gthth = 1./x[1]**2
        gphph = 1./(x[1]*sin(x[2]))**2
        gtph = 0.
        
        return [gtt, grr, gthth, gphph, gtph]
    
    def dr_inverse_metric(self,x):
        '''
        This procedure returns the derivative of the Schwarzschild inverse metric 
        w.r.t the coordinate r
        ===========================================================================
        Coordinates 
        t = x[0]
        r = x[1]
        theta = x[2]
        phi = x[3]
        ===========================================================================
        '''
        # Derivative of the metric components
        drgtt =  self.drgtt(x[1])
        drgrr = self.drgrr(x[1])
        drgthth = -2/x[1]**3
        drgphph = -2/(x[1]**3*sin(x[2])**2)
        drgtph = 0.
        
        return [drgtt, drgrr, drgthth, drgphph, drgtph]


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
        # Metric and its numerical derivative
        gtt, grr, gthth, gphph, gtph = self.inverse_metric(q[0:4])
        drgtt, drgrr, drgthth, drgphph, drgtph = self.dr_inverse_metric(q[0:4])
        
        # Geodesics differential equations 
        dtdlmbda = gtt*q[4]
        drdlmbda = grr*q[5]
        dthdlmbda = gthth*q[6]
        dphidlmbda = gphph*q[7]
        
        dk_tdlmbda = 0.
        dk_rdlmbda = - (drgtt*q[4]**2)/2 - (drgrr*q[5]**2)/2 \
                     - (drgthth*q[6]**2)/2 - (drgphph*q[7]**2)/2
        dk_thdlmbda = (cos(q[2])/sin(q[2])**3)*(q[7]/q[1])**2
        dk_phidlmbda = 0.
        
        return [dtdlmbda, drdlmbda, dthdlmbda, dphidlmbda, 
                dk_tdlmbda, dk_rdlmbda, dk_thdlmbda, dk_phidlmbda]





###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')

    

