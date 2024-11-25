# PLAN
# Construct Sparse Matrix either using function above or similar way
# Solve eigenvalues/vectors and check against inspo
# Plot
# Optimise
# Parallise

# Optional Features: adding the boundary back on

import argparse
import sys
import numpy as np
from PIL import Image
import scipy as sp
from Debugging import Debugging


# Function to get potential from image
def get_potential_from_image(file_name: str, debugger: Debugging) -> tuple:

    # Open image in black and white
    image = (Image.open(file_name)).convert('1')

    # Create an image matrix with white pixels being false, and black being true
    bool_image_matrix = np.array(image) == 0

    # Get number of points in y and x (num of rows/cols in the original image)
    total_y_points, total_x_points = bool_image_matrix.shape

    # Image needs to be a square to construct a square Hamiltonian
    if total_x_points != total_y_points:
        sys.exit("ERROR: Image must be a square")

    # Create a mesh for when constructing the Hamiltonian
    potential_matrix = np.where(bool_image_matrix, 0, 1)

    # Save debug information
    debugger.debug_store(bool_image_matrix, "./debug/bool_image.mat")
    debugger.debug_store(potential_matrix, "./debug/position_mesh.mat")
    
    return (potential_matrix.flatten(), (total_x_points, total_y_points))

def construct_sparse_hamiltonian(potential_info: tuple, debugger: Debugging) -> sp.sparse.csr:
    """Constructs a sparse Hamiltonian matrix."""
    potential_matrix = potential_info[0] # Get relevant information 
    num_of_points = potential_info[1][0] # This is total of x-points
    N = num_of_points**2 # Hamiltonian needs NxN as each physical point needs a row in overall matrix
    x_intervals = np.linspace(0, 1, num_of_points) # Change this to work on any grid size
    grid_step = np.abs(x_intervals[1] - x_intervals[0])
    # Need for batch population
    rows, cols, data = [], [], []

    denom = grid_step**2 # Grid step squared for Forward Difference
    denom_inv = -1 / denom # Avoids repeated calculation

    # Only need to focus on the terms that will be non-zero
    for j in range(N):
        # Only consider the Hamiltonian where the potential is 0
        if potential_matrix[j] == 0:

            # Diaganol term will always be non-zero
            rows.append(j)
            cols.append(j)
            data.append(-4*denom_inv)

            # Get neighbour info
            right = j + 1
            left = j - 1
            up = j + num_of_points
            down = j - num_of_points

            # This is the Forward Diff method
            # Checks the nearest neighbours and if they are within the potential and not on the boundaries
            if  right < N and right % num_of_points != 0 and potential_matrix[right] == 0:
                rows.append(j)
                cols.append(right)
                data.append(denom_inv)

            if left >= 0 and j % num_of_points != 0 and potential_matrix[left] == 0:
                rows.append(j)
                cols.append(left)
                data.append(denom_inv)

            if down >= 0 and potential_matrix[down] == 0:
                rows.append(j)
                cols.append(down)
                data.append(denom_inv)

            if up < N and potential_matrix[up] == 0:
                rows.append(j)
                cols.append(up)
                data.append(denom_inv)

    # Returns Hamiltonian by creating sparse matrix in batch
    return sp.sparse.coo_matrix((data, (rows, cols)), shape=(N, N)).tocsr() # Return compressed sparse row matrix for eigen solver
