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