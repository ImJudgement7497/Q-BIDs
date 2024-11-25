# Use sparse matrices instead of dense matrices 
# Basically exact code for getting the image, just made more robust.
# Try not to use Finite Difference, and implement other method

# from scipy.sparse import lil_matrix

# def construct_sparse_hamiltonian(N, increment, position_mesh):
#     """Constructs a sparse Hamiltonian matrix."""
#     H = lil_matrix((N, N))

#     for j in range(N):
#         if position_mesh[j] != 0:
#             denom = increment**2
#             H[j, j] = 4 / denom

#             if j + 1 < N and (j + 1) % temp != 0 and position_mesh[j + 1] != 0:
#                 H[j, j + 1] = -1 / denom

#             if j - 1 >= 0 and j % temp != 0 and position_mesh[j - 1] != 0:
#                 H[j, j - 1] = -1 / denom

#             if j - temp >= 0 and position_mesh[j - temp] != 0:
#                 H[j, j - temp] = -1 / denom

#             if j + temp < N and position_mesh[j + temp] != 0:
#                 H[j, j + temp] = -1 / denom

#     return H

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
from Debugging import Debugging
import scipy as sp

# Function to get potential from image
def get_potential_from_image(file_name: str) -> tuple: 

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

def construct_sparse_hamiltonian(N: int, increment: float, potential_matrix: np.ndarray) -> sp.sparse.csr:
    """Constructs a sparse Hamiltonian matrix."""
    H = sp.sparse.lil_matrix((N, N))

    # Only need to focus on the terms that will be non-zero
    for j in range(N):
        # Only consider the Hamiltonian where the potential is 0
        if potential_matrix[j] == 0:

            denom = increment**2
            temp = int(np.sqrt(N))

            # Diaganol term will always be non-zero
            H[j, j] = 4 / denom

            if j + 1 < N and (j + 1) % temp != 0 and potential_matrix[j + 1] == 0:
                H[j, j + 1] = -1 / denom

            if j - 1 >= 0 and j % temp != 0 and potential_matrix[j - 1] == 0:
                H[j, j - 1] = -1 / denom

            if j - temp >= 0 and potential_matrix[j - temp] == 0:
                H[j, j - temp] = -1 / denom

            if j + temp < N and potential_matrix[j + temp] == 0:
                H[j, j + temp] = -1 / denom

    return H.tocsr()

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true', help='Enable debug mode')
args = parser.parse_args()

debugger = Debugging(debug=args.debug)

temp = get_potential_from_image(file_name="./three.png")
Hamiltonian = construct_sparse_hamiltonian(temp[1][0] * temp[1][1], 1, temp[0])

eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=8)

