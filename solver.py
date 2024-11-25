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

#             if j + 1 < N and (j + 1) % int(np.sqrt(N)) != 0 and position_mesh[j + 1] != 0:
#                 H[j, j + 1] = -1 / denom

#             if j - 1 >= 0 and j % int(np.sqrt(N)) != 0 and position_mesh[j - 1] != 0:
#                 H[j, j - 1] = -1 / denom

#             if j - int(np.sqrt(N)) >= 0 and position_mesh[j - int(np.sqrt(N))] != 0:
#                 H[j, j - int(np.sqrt(N))] = -1 / denom

#             if j + int(np.sqrt(N)) < N and position_mesh[j + int(np.sqrt(N))] != 0:
#                 H[j, j + int(np.sqrt(N))] = -1 / denom

#     return H

# PLAN
# Construct Sparse Matrix either using function above or similar way
# Solve eigenvalues/vectors and check against inspo
# Plot
# Optimise
# Parallise

# Optional Features: adding the boundary back on

import argparse
import os
import pickle
import numpy as np
from PIL import Image
from Debugging import Debugging

# Function to get potential from image
def get_potential_from_image(file_name: str) -> tuple: 

    # Open image in black and white
    image = (Image.open(file_name)).convert('1')

    # Create an image matrix with white pixels being false, and black being true
    bool_image_matrix = np.array(image) == 0

    # Get number of points in y and x (num of rows/cols in the original image)
    total_y_points, total_x_points = bool_image_matrix.shape

    # Create a mesh for when constructing the Hamiltonian
    position_mesh_matrix = np.where(bool_image_matrix, 0, 1)

    # Save debug information
    debugger.debug_store(bool_image_matrix, "./debug/bool_image.mat")
    debugger.debug_store(position_mesh_matrix, "./debug/position_mesh.mat")
    
    return (position_mesh_matrix.flatten(), (total_x_points, total_y_points))

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true', help='Enable debug mode')
args = parser.parse_args()

debugger = Debugging(debug=args.debug)

temp = get_potential_from_image(file_name="./three.png")