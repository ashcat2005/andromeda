"""
===============================================================================
Main script 
Creates the Black Hole image
===============================================================================
@author: Eduard Larrañga - 2023
===============================================================================
"""
#import warnings
#warnings.filterwarnings('ignore')

from numpy import pi
#from metrics import schwarzschild
from metrics import kerr
from accretion_structures import thin_disk 
from detectors import image_plane 
from common.common import Image





###############################################################################
############################ BLACK HOLE DEFINITION ############################
###############################################################################
M = 1 # Mass
a = 0.5 # Angular Monmentum
#blackhole = schwarzschild.BlackHole(M)
blackhole = kerr.BlackHole(M,a)


###############################################################################
############################# DETECTOR PARAMETERS #############################
###############################################################################
D = 100*M
iota = pi/2.3
x_screen_side = 25*M
y_screen_side = 15*M
n_pixels = 5
detector = image_plane.detector(D=D, iota = iota, 
                       x_s_side = x_screen_side, y_s_side = y_screen_side,
                       n_pixels=n_pixels)


###############################################################################
############################# ACCRETION STRUCTURE #############################
###############################################################################
#R_min = blackhole.ISCOco 
#R_max = 20*M
acc_structure = thin_disk.structure(blackhole)#, R_min = 10)


###############################################################################
############################### IMAGE FILENAME ################################
###############################################################################
filename = 'NovikovThorneBlackHole'
savefig = False









#################################### MAIN #####################################

image = Image()

# Photons creation
image.create_photons(blackhole, detector)

# Create the image data
image.create_image(blackhole, acc_structure)

# Plot the image
image.plot(savefig=savefig, filename=filename)
