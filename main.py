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
from metrics import schwarzschild
from metrics import kerr
from metrics import num_schwarzschild
from metrics import scalar_hair_BH
from accretion_structures import simple_disk 
from accretion_structures import thin_disk 
from accretion_structures import ring_disk 
from accretion_structures import ring_disk_2
from detectors import image_plane 
from common.common import Image




'''
===============================================================================
============================ BLACK HOLE DEFINITION ============================
===============================================================================
'''
##### SCHWARZSCHILD BH
#M = 1 # Mass
#blackhole = schwarzschild.BlackHole(M)

##### KERR BH
M = 1 # Mass
a = 0.4 # Angular Monmentum
blackhole = kerr.BlackHole(M,a)

##### NUMERICAL SCHWARZSCHILD BH
#blackhole = num_schwarzschild.BlackHole(M)

##### SCALAR HAIR BH
#M = 1
#blackhole = scalar_hair_BH.BlackHole()

'''
===============================================================================
=========================== DETECTOR PARAMETERS ===============================
===============================================================================
'''
D = 100*M
iota = (85)*pi/180 #pi/2.3
x_screen_side = 25*M
y_screen_side = 20*M
n_pixels = 2
detector = image_plane.detector(D=D, iota = iota, 
                       x_s_side = x_screen_side, y_s_side = y_screen_side,
                       n_pixels=n_pixels)


'''
===============================================================================
============================ ACCRETION STRUCTURE ==============================
===============================================================================
'''

#acc_structure = simple_disk.structure(blackhole, 6, 20)


########### NOVIKOV-THORNE THIN DISK
R_min = blackhole.ISCOco 
R_max = 20*M
acc_structure = thin_disk.structure(blackhole, R_min = 6)

########### RINGS DISK

#acc_structure = ring_disk.structure(blackhole)


'''
===============================================================================
============================ IMAGE FILENAME ===================================
===============================================================================
'''
filename = 'ScalarHairBlackHoleNT3'
savefig = False








'''
===============================================================================
==================================== MAIN =====================================
===============================================================================
'''
image = Image()

# Photons creation
image.create_photons(blackhole, detector)

# Create the image data
image.create_image(blackhole, acc_structure)

# Plot the image
image.plot(savefig=savefig, filename=filename, cmap='inferno')
#image.plotContours(savefig=savefig, filename=filename)
