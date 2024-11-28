import numpy as np
import pickle
import scipy as sp

with open("./debug/Ham.mat.pkl", "rb") as file:
    Hamiltonian = pickle.load(file)
    Hamiltonian = Hamiltonian.toarray()

for row in Hamiltonian:
    print(row)