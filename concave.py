import numpy as np
import matplotlib.pyplot as plt

# Define the implicit function
def implicit_function(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return (r - 2)**2 + (theta**2) / 4 - 1  # The equation we derived

# Define the grid range
x_vals = np.linspace(-4, 4, 400)
y_vals = np.linspace(-4, 4, 400)
X, Y = np.meshgrid(x_vals, y_vals)

# Compute function values
Z = implicit_function(X, Y)

# Define the potential matrix
V0 = 1e6  # Example: Set outside potential to 10
V = np.where(Z <= 0, 0, V0)  # 0 inside, V0 outside

# Plot the potential matrix
plt.figure(figsize=(6, 6))
plt.imshow(V, extent=[-4, 4, -4, 4], origin='lower', cmap='coolwarm', alpha=0.8)
plt.colorbar(label="Potential V(x, y)")
plt.contour(X, Y, Z, levels=[0], colors='black')  # Boundary outline
plt.xlabel("x")
plt.ylabel("y")
plt.title("Potential Matrix")
plt.grid()

# Save and show the plot
plt.savefig("potential_matrix.png", dpi=300)
plt.show()
