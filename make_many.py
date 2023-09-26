"""
===============================================================================
This script produces many images of Kerr black holes with different values of 
the spin parameter and/or inclination angles.
===============================================================================
@author: Eduard Larrañga - 2023
===============================================================================
"""
#import warnings
#warnings.filterwarnings('ignore')

from numpy import pi
from black_holes import schwarzschild
from black_holes import kerr
from black_holes import num_schwarzschild
from black_holes import scalar_hair_BH
from accretion_structures import simple_disk 
from accretion_structures import thin_disk 
from accretion_structures import ring_disk 
from accretion_structures import ring_disk_2
from detectors import image_plane 
from common.common import Image
import numpy as np



############# Common Parameters

##### KERR BH
M = 1   # Mass

D = 100*M  # Distance to the BH
x_side = 25*M
x_pixels = 90
y_pixels = int(x_pixels*9/16)


N = 1000                    # Number of images
aa = np.random.random(N)    # Angular Monmentum

images_data = np.zeros([N,x_pixels,y_pixels])
labels = np.zeros([N])


for i in range(N):
    print('\nImage number: ', i+1)
    a = aa[i]
    blackhole = kerr.BlackHole(M,a)

    iota = (pi/180)*(80)   # Inclination Angle
    detector = image_plane.detector(D=D, iota=iota, x_pixels=x_pixels, 
                                    x_side=x_side, ratio='16:9')

    #acc_structure = simple_disk.structure(blackhole, 6, 20)


    ########### NOVIKOV-THORNE THIN DISK
    #R_min = blackhole.ISCOco 
    #R_max = 20*M
    acc_structure = thin_disk.structure(blackhole)

    image = Image()
    # Photons creation
    image.create_photons(blackhole, detector)
    # Create the image data
    image.create_image(blackhole, acc_structure)

    labels[i] = a
    images_data[i] = image.image_data

    if i%50==0:
        print('\nSaving info ...')
        np.save('images_data.npy', images_data)
        np.save('labels.npy', labels)


np.save('images_data.npy', images_data)
np.save('labels.npy', labels)
