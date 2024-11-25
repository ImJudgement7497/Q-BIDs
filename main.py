#!/usr/bin/env python3

import scipy as sp
from Debugging import Debugging
from Parser import Parser
from solver import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import os

parser = Parser("./input.txt")

image_file = parser.get_config_value("potential_image_name")
if not os.path.isdir(f"./{image_file}_plots"):
    os.mkdir(f"./{image_file}_plots")
debugger = Debugging(debug=parser.args.debug)

potential_info = get_potential_from_image(file_name=image_file, debugger=debugger)
Hamiltonian = construct_sparse_hamiltonian(potential_info, debugger=debugger)

# Eigenvalues and eigenvectors from eigsh
eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=1)

idx = eigenvalues.argsort()  # Sort in ascending order
print(eigenvalues[idx])
eigenvalues = eigenvalues[idx]  # Sorted eigenvalues
idx = eigenvalues.argsort()[::-1]
eigenvectors = eigenvectors[:, idx]  # Reordered eigenvectors to match eigenvalues

print("Eigenvalues (sorted):", eigenvalues)

# Define M, N based on Hamiltonian size
M = int(np.sqrt(Hamiltonian.shape[0]))  # Assuming Hamiltonian is square
N = M  # Adjust if the system is non-square

# Loop through eigenvectors (eigenvectors are returned in columns)
for i in range(eigenvectors.shape[1]):  # iterate over each eigenvector (column)
    figi, axi = plt.subplots(1, 1)
    
    # Reshape eigenvector column 'i' to MxN grid and calculate its squared absolute value
    eigenvector_reshaped = np.transpose(np.absolute(eigenvectors[:, i].reshape(M, N)) ** 2)
    
    # Plot the reshaped eigenvector
    plot = plt.imshow(np.transpose(eigenvector_reshaped), cmap='magma', interpolation='gaussian')
    
    # Remove axis ticks
    plt.setp(axi, xticks=[], yticks=[])
    
    # Add colorbar
    divider = make_axes_locatable(axi)
    cax = divider.append_axes("right", size="3%", pad=0.1)
    cbar = figi.colorbar(plot, ax=axi, extend='both', cax=cax)
    cbar.minorticks_on()
    cbar.ax.tick_params(labelsize=5, pad=0.1)
    
    # Set title based on the eigenstate
    if i == 0:
        axi.set_title('The ground state', fontsize=12)
    elif i == 1:
        axi.set_title('The 1$^{st}$ excited state', fontsize=12)
    elif i == 2:
        axi.set_title('The 2$^{nd}$ excited state', fontsize=12)
    elif i == 3:
        axi.set_title('The 3$^{rd}$ excited state', fontsize=12)
    else:
        axi.set_title(f'{i}$^{{th}}$ excited state', fontsize=12)
    
    # Label the eigenvalue on the plot
    axi.text(0.5, 0.95, f'Eigenvalue: {eigenvalues[i]:.3f}', transform=axi.transAxes, ha='center', fontsize=10, color='white')
    
    # Save the plot to a PDF file
    plt.savefig(f'./{image_file}_plots/{i}.png')
    plt.close(figi)


