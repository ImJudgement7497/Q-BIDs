import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
from scipy import sparse
from numba import jit
import os
my_path = os.path.abspath(r"./paper")
##Create meshgrid for x and y
N = 100
L = 1
X,Y = np.meshgrid(np.linspace(0,L,N,
dtype=float),np.linspace(0,L,N,
dtype=float))
##potential m\deltax^2 unit
p = "2DHA"
V = 0*X
##create matrix
diag = np.ones([N])

diags = np.array([diag, -2*diag, diag])
D = sparse.spdiags(diags, np.array([-1, 0,1]), N, N)
##define energy
T = -1/2 * sparse.kronsum(D,D)
U = sparse.diags(V.reshape(N**2),(0))
H = T+U
##Solve for eigenvector and eigenvalue
eigenvalues , eigenvectors = eigsh(H, k=10, which='SM')
# for row in H.toarray():
#     print(row)
print(eigenvalues)

def get_e(n):
    return eigenvectors.T[n].reshape((N,N))
##number of state
# for n in range (0,10):
#     ##plot V
#     plot0 = plt.figure(0,figsize=(8,6))
#     cs = plt.contourf(X,Y,V,100)
#     plt.colorbar()
#     for c in cs.collections:
#         c.set_rasterized(True)
#     plt.title("Plot of V")
#     plt.xlabel(r"$X$")
#     plt.ylabel(r"$Y")
#     plt.savefig(os.path.join(my_path,
#     "Figure_{}.{}.0.pdf".format(p,n)))
#     ##plot eigenvector
#     plot1 = plt.figure(1,figsize=(8,6))
#     cs = plt.contourf(X,Y, get_e(n),300)
#     plt.colorbar()
#     for c in cs.collections:
#         c.set_rasterized(True)
#     plt.title("Plot of Eigenfunction for {} state".format(n))
#     # plt.xlabel(r"$X$")
#     # plt.ylabel(r"$Y$"s)
#     plt.savefig(os.path.join(my_path,
#     "Figure_{}.{}.1.pdf".format(p,n)))
#     ##plot probability density
#     plot2 = plt.figure(2, figsize=(8,6))
#     cs = plt.contourf(X,Y, get_e(n)**2,300)
#     plt.colorbar()
#     for c in cs.collections:
#         c.set_rasterized(True)
#     plt.title("Plot of Probability Density for {} state".format(n))
#     # plt.xlabel(r’$X$’)
#     # plt.ylabel(r’$Y$’)
#     plt.savefig(os.path.join(my_path,
#     "Figure_{}.{}.2.pdf".format(p,n)))
##plot eigenvalues
plot3 = plt.figure(3)
alpha = eigenvalues[0]/2
E_a = eigenvalues /alpha
b = np.arange(0, len(eigenvalues),1)
plt.scatter(b, E_a, s=1444, marker="_", linewidth=2, zorder=3)
plt.title("Plot of eigenvalues")
plt.savefig("./paper/eigenvalues.png")
# plt.xlabel(’$(n_{x})^2+(n_{y})^2$’)
# plt.ylabel(r’$mE/\hbar^2$’)