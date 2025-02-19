#!/usr/bin/env python3

import scipy as sp
from Debugging import Debugging
from Parser import Parser
from solver import *
from plotting import *
import function
import matplotlib.pyplot as plt
import os

parser = Parser("./input.txt")

if not os.path.isdir("./plots"):
    os.mkdir("./plots")
shape = parser.get_config_value("shape")
image_file = parser.get_config_value("potential_image_name")
function_parsed = parser.get_config_value("function")

if shape != None and image_file != None:
    sys.exit("Please provide only a shape or a potential_image_name")
if shape == None and image_file == None and function_parsed == None:
    sys.exit("Please provide a shape, a potential_image_name or a function")

function_bool = False

if function_parsed != None:
    function_bool = True
    
if shape != None and function_bool == False:
    potential_type = shape
    is_shape = True
    print(f"Using shape {potential_type}")
elif shape == None and function_bool == False:
    potential_type = image_file
    is_shape = False
    print(f"Using file {potential_type}")
else:
    potential_type = function_parsed
    is_shape = True # Needed for plotting, I know it does not make sense
    print("Using function defined in function.py")

if not os.path.isdir(f"./plots/{potential_type}_plots"):
    os.mkdir(f"./plots/{potential_type}_plots")

if not os.path.isdir("./results"):
    os.mkdir("./results")

max_level = parser.get_config_value("max_level")
if max_level == None:
    sys.exit("Please provide max_levels")

for i in range(max_level):
    if not os.path.isdir(f"./plots/{potential_type}_plots/{i}"):
        os.mkdir(f"./plots/{potential_type}_plots/{i}")

boundary_value = parser.get_config_value("boundary_value")
if boundary_value == None:
    boundary_value = 1e6

grid_info_names = ["x_start", "x_end", "y_start", "y_end", "grid_size"]
grid_info = []
for elm in grid_info_names:
    temp = parser.get_config_value(elm)
    grid_info.append(temp)

grid_step = 1 / (grid_info[4]-1)
grid_info.append(grid_step)

debugger = Debugging(debug=parser.args.debug)

if function_bool:
    x_vals = np.linspace(grid_info[0], grid_info[1], grid_info[4])
    y_vals = np.linspace(grid_info[0], grid_info[1], grid_info[4])
    X, Y = np.meshgrid(x_vals, y_vals)
    
    Z = function.function(X, Y)
    
    plt.figure(figsize=(8, 6))
    plt.contour(X, Y, Z, levels=[0], colors='blue')
    plt.title("Function in Cartesian Coordinates")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig(f"./plots/{potential_type}_plots/function_plot.png")
    potential_matrix = np.where(Z <= 0, 0, boundary_value)
    
    if potential_matrix.size != grid_info[4] * grid_info[4]:
        raise ValueError(f"Expected {grid_info[4]*grid_info[4]} elements, but got {potential_matrix.size}.")
    potential_matrix = potential_matrix.reshape(grid_info[4], grid_info[4])

    
    Hamiltonian = construct_hamiltonian(potential_matrix, grid_info[4], debugger, boundary_value)
    print("Solving now")
    eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=max_level, which="SM")
    eigenvectors = eigenvectors * (grid_step**2)
    plot_potential(potential_matrix, potential_type, grid_info, is_shape)
    plot_eigenfunctions_from_shape(eigenvectors, grid_info[4], max_level, potential_type)
    plot_nodal_lines(eigenvectors, grid_info[4], max_level, potential_type)
    
else:
    
    if not is_shape:
        potential_info = get_potential_from_image(image_file, debugger, boundary_value)
        Hamiltonian = construct_hamiltonian(potential_info[0], potential_info[1], debugger, boundary_value)
        print("Solving now")
        eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=max_level, which="SM")
        eigenvectors = eigenvectors * (grid_step**2)
        plot_potential(potential_info[0], potential_type, grid_info, is_shape)
        plot_eigenfunctions_from_image(eigenvectors, potential_info, max_level, potential_type)
        plot_nodal_lines(eigenvectors, potential_info[1], max_level, potential_type)
    else:
        potential_matrix = get_potential_from_shape(shape, grid_info, debugger, boundary_value)
        Hamiltonian = construct_hamiltonian(potential_matrix, grid_info[4], debugger, boundary_value)
        print("Solving now")
        eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=max_level, which="SM")
        eigenvectors = eigenvectors * (grid_step**2)
        plot_potential(potential_matrix, potential_type, grid_info, is_shape)
        plot_eigenfunctions_from_shape(eigenvectors, grid_info[4], max_level, potential_type)
        plot_nodal_lines(eigenvectors, grid_info[4], max_level, potential_type)


# Eigenvalues and eigenvectors from eigsh

np.save(f"./results/{potential_type}_eigenvectors_upto_state_{max_level}.npy", eigenvectors)
np.save(f"./results/{potential_type}_eigenvalues_upto_state_{max_level}.npy", eigenvalues)

# plot_eigenfunctions(eigenvectors, potential_info, max_level, image_file)
# plot_prob_densities(eigenvectors, potential_info, max_level, image_file)
# plot_eigenfunction_zero_crossings(eigenvectors, potential_info, max_level, image_file)
