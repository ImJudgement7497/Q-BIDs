#!/usr/bin/env python3

import scipy as sp
from Debugging import Debugging
from Parser import Parser
from solver import *
from plotting import *
import matplotlib.pyplot as plt
import os

parser = Parser("./input.txt")

if not os.path.isdir("./plots"):
    os.mkdir("./plots")

image_file = parser.get_config_value("potential_image_name")
if not os.path.isdir(f"./plots/{image_file}_plots"):
    os.mkdir(f"./plots/{image_file}_plots")

if not os.path.isdir("./results"):
    os.mkdir("./results")

max_level = parser.get_config_value("max_level")
if max_level == None:
    sys.exit("Please provide max_levels")

for i in range(max_level):
    if not os.path.isdir(f"./plots/{image_file}_plots/{i}"):
        os.mkdir(f"./plots/{image_file}_plots/{i}")

boundary_value = parser.get_config_value("boundary_value")
if boundary_value == None:
    boundary_value = 1e6

debugger = Debugging(debug=parser.args.debug)

potential_info = get_potential_from_image(image_file, debugger, boundary_value)
Hamiltonian = construct_hamiltonian(potential_info, debugger, boundary_value)

N = potential_info[1]
grid_step = 1 / (N-1)
eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=max_level, which="SM")
eigenvectors = eigenvectors * (grid_step**2)
# Eigenvalues and eigenvectors from eigsh

np.save(f"./results/{image_file}_eigenvectors_upto_state_{max_level}.npy", eigenvectors)
np.save(f"./results/{image_file}_eigenvalues_upto_state_{max_level}.npy", eigenvalues)

plot_potential(potential_info, image_file)
plot_eigenfunctions(eigenvectors, potential_info, max_level, image_file)
plot_prob_densities(eigenvectors, potential_info, max_level, image_file)
# plot_eigenfunction_zero_crossings(eigenvectors, potential_info, max_level, image_file)
