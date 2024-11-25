import scipy as sp
from Debugging import Debugging
from Parser import Parser
from solver import *

parser = Parser("./input.txt")

image_file = parser.get_config_value("potential_image_name")
debugger = Debugging(debug=parser.args.debug)

potential_info = get_potential_from_image(file_name=image_file, debugger=debugger)
Hamiltonian = construct_sparse_hamiltonian(potential_info, debugger=debugger)

eigenvalues, eigenvectors = sp.sparse.linalg.eigsh(Hamiltonian, k=8)