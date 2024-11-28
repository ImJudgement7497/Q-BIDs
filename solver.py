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
    
    N = total_x_points

    # Create a mesh for when constructing the Hamiltonian

    potential_matrix = np.where(bool_image_matrix, 0, 1) # MAKE THIS USER INPUTTED

    # Save debug information
    debugger.debug_store(bool_image_matrix, "./debug/bool_image.mat")
    debugger.debug_store(potential_matrix, "./debug/position_mesh.mat")
    
    return (potential_matrix, N)

def construct_hamiltonian(potential_info, debugger: Debugging, boundary_value=1000):

    potential_matrix = potential_info[0]
    N = potential_info[1]

    r"""
    Following created by Wai Jui Wong, paper at https://www.researchgate.net/publication/356858518_Solving_2D_Time_Independent_Schrodinger_Equation_Using_Numerical_Method
    """
    #Construct 2D Laplacian
    diag = np.ones([N])
    diags = np.array([diag, -2*diag, diag])
    D = sp.sparse.spdiags(diags, np.array([-1, 0, 1]), N, N)

    # Construct Kinetic energy operator
    # See paper for why we do this
    K = -1/2 * sp.sparse.kronsum(D, D)

    V = potential_matrix.reshape(N, N)  # Reshape into 2D grid
    
    # Apply infinite potential at the boundaries by setting large values at the boundary indices
    V[0, :] = boundary_value  # Top boundary
    V[-1, :] = boundary_value  # Bottom boundary
    V[:, 0] = boundary_value  # Left boundary
    V[:, -1] = boundary_value  # Right boundary
    
    # Convert the potential into a sparse matrix
    V_sparse = sp.sparse.diags(V.flatten(), 0)  # Flatten V and convert to sparse

    # print(V_sparse.toarray())

    # Construct H
    H = K + V_sparse

    debugger.debug_store(H, "./debug/Ham.mat")
    # for row in H.toarray():
    #     print(row)
    return H