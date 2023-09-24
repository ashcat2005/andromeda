from mpi4py import MPI

def init_mpi():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    return comm, rank, size


def distribute(rank, size, detector):
    total_photons = detector.x_pixels * detector.y_pixels
    photons_per_rank = total_photons // size 
    rest = total_photons % size

    if rank < rest:
        start_idx = rank * (photons_per_rank + 1)
        end_idx = start_idx + photons_per_rank + 1
    else:
        start_idx = (rank * photons_per_rank) + rest
        end_idx = start_idx + photons_per_rank

    start_alpha, start_beta = divmod(start_idx, len(detector.betaRange))
    end_alpha, end_beta = divmod(end_idx, len(detector.betaRange))

    return start_alpha, start_beta, end_alpha, end_beta