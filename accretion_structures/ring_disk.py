"""
===============================================================================
Model of concentrical rings for the accretion disk 
===============================================================================
@author: Eduard Larrañaga - 2023
===============================================================================
"""

class structure:
    def __init__(self, blackhole, corotating=True):
        self.in_edge = 5.9
        self.out_edge = 20

    def energy_flux(self, r):
        tol = 0.2
        radii = [6,8,10,12,14,16,18]
        intensity = 0.
        for radius in radii:
            if r>radius-tol and r<radius+tol:
                intensity = 1.
        return intensity



###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')
