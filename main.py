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

from numpy import pi, zeros, sum
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

from parallel import *


'''
===============================================================================
============================ BLACK HOLE DEFINITION ============================
===============================================================================
'''

##### SCHWARZSCHILD BH
#M = 1 # Mass
#blackhole = schwarzschild.BlackHole(M)


##### KERR BH
M = 1   # Mass
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
D = 100*M              # Distance to the BH
iota = (pi/180)*(85)   # Inclination Angle
x_screen_side = 25*M
y_screen_side = 20*M
n_pixels = 5
detector = image_plane.detector(D=D, iota = iota, 
                                x_s_side = x_screen_side, 
                                y_s_side = y_screen_side,
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
filename = 'Parallel-ver'
savefig = False


'''
===============================================================================
==================================== MAIN =====================================
===============================================================================
'''
# Measure time
start_time = MPI.Wtime()

#initialize parallelization variables
comm, rank, size = init_mpi()
start_alpha, start_beta, end_alpha, end_beta = distribute(rank, size, detector)


image = Image(start_alpha, start_beta, end_alpha, end_beta)

# Photons creation
image.create_photons(blackhole, detector, rank)

# Create the image data
image.create_image(blackhole, acc_structure,rank)


gathered_image_data = comm.gather(image.image_data, root=0)

# Check if the current process is the principal
if rank == 0:
    
    final_image_data = zeros([detector.x_pixels, detector.y_pixels])

    # Sum up the image data gathered from all processes to get the final image.
    final_image_data = sum(gathered_image_data, axis=0)

    # Assign the final image data to the respective variable.
    image.image_data = final_image_data

    end_time = MPI.Wtime()
    
    #total elapsed time.
    elapsed_time = end_time - start_time
    print(f"\n TIME: {elapsed_time} sec")

    # Plot the image using the final image data.
    image.plot(savefig=savefig, filename=filename, cmap='inferno')
    # image.plotContours(savefig=savefig, filename=filename)
