"""
===============================================================================
Parameters to create the Black Hole image
===============================================================================
@author: Eduard Larrañga - 2023
===============================================================================
"""

from numpy import pi
#from metrics import schwarzschild
from metrics import kerr
from accretion_structures import thin_disk 
from image_plane import *



###############################################################################
############################ BLACK HOLE DEFINITION ############################
###############################################################################
M = 1 # Mass
a = 0.8 # Angular Monmentum
#blackhole = schwarzschild.BlackHole(M)
blackhole = kerr.BlackHole(M,a)


###############################################################################
############################# DETECTOR PARAMETERS #############################
###############################################################################
D = 100*M
iota = pi/2.1
x_screen_side = 25*M
y_screen_side = 15*M
n_pixels = 10
detector = image_plane(D=D, iota = iota, 
                       x_s_side = x_screen_side, y_s_side = y_screen_side,
                       n_pixels=n_pixels)


###############################################################################
############################# ACCRETION STRUCTURE #############################
###############################################################################
R_min = blackhole.ISCOco 
R_max = 20*M
acc_structure = thin_disk.structure(blackhole, R_min, R_max)


###############################################################################
############################### IMAGE FILENAME ################################
###############################################################################
filename = 'RotatingBlackHoleThinDisk.png'
savefig = True



###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')
