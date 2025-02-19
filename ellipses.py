import numpy as np
import matplotlib.pyplot as plt

# # Define the implicit function (adjustable for multiple contours)
# def implicit_function(x, y):
#     c = 4
#     B = x**2 + y**2 - c**2
#     p = (-B + np.sqrt((B**2) + 4 * (c**2) * (y**2))) / (2 * (c**2))
#     q = (-B - np.sqrt((B**2) + 4 * (c**2) * (y**2))) / (2 * (c**2))
    
#     eta_0 = np.arcsin(np.sqrt(p))
    
#     # Initialize eta with zeros
#     eta = np.zeros_like(x)
    
#     # Quadrant I (x >= 0, y >= 0)
#     eta[(x >= 0) & (y >= 0)] = eta_0[(x >= 0) & (y >= 0)]
    
#     # Quadrant II (x < 0, y >= 0)
#     eta[(x < 0) & (y >= 0)] = np.pi - eta_0[(x < 0) & (y >= 0)]
    
#     # Quadrant III (x =< 0, y < 0)
#     eta[(x <= 0) & (y < 0)] = np.pi + eta_0[(x <= 0) & (y < 0)]
    
#     # Quadrant IV (x > 0, y < 0)
#     eta[(x > 0) & (y < 0)] = 2 * np.pi - eta_0[(x > 0) & (y < 0)]
    
#     xi = 0.5 * np.log(1 - 2*q + 2 * np.sqrt((q**2) - q))
    
#     return xi**2 + eta**2 - 8

# # Define grid
# x_vals = np.linspace(-10, 10, 400)
# y_vals = np.linspace(-10, 10, 400)
# X, Y = np.meshgrid(x_vals, y_vals)

# Z = implicit_function(X, Y)

# # Define the potential matrix
# V0 = 10  # Example: Set outside potential to 10
# V = np.where(Z <= 0, 0, V0)  # 0 inside, V0 outside

# # Plot the potential matrix
# plt.figure(figsize=(6, 6))
# # plt.imshow(V, extent=[-4, 4, -4, 4], origin='lower', cmap='coolwarm', alpha=0.8)
# # plt.colorbar(label="Potential V(x, y)")
# plt.contour(X, Y, Z, levels=[0], colors='black')  # Boundary outline
# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("Potential Matrix")
# plt.grid()

# # Save and show the plot
# plt.savefig("ELLIPTICAL.png", dpi=300)
# plt.show()

# # Function for the superellipse in Cartesian coordinates
# def superellipse(x, y, a, b, n):
#     return (np.abs(x) / a)**n + (np.abs(y) / b)**n - 1

# # Generate a grid for x and y
# x_vals = np.linspace(-5, 5, 400)
# y_vals = np.linspace(-5, 5, 400)
# X, Y = np.meshgrid(x_vals, y_vals)

# # Parameters for the superellipse
# a = 1
# b = 2
# n = 1.5 # Change this value to adjust the convexity/concavity

# # Evaluate the superellipse
# Z = superellipse(X, Y, a, b, n)

# # Plot the superellipse in Cartesian coordinates (contour)
# plt.figure(figsize=(8, 6))
# contour = plt.contour(X, Y, Z, levels=[0], colors='blue')
# plt.title("Superellipse in Cartesian Coordinates")
# plt.xlabel("X")
# plt.ylabel("Y")
# plt.show()

# Superellipse equation in elliptical coordinates
def superellipse_elliptical(eta, a, b, n):
    term1 = (np.cosh(eta) / a) ** n
    term2 = (np.sinh(eta) / b) ** n
    xi = (1 / (term1 + term2)) ** (1 / n)
    return xi

# Parameters for the superellipse
a = 4  # semi-major axis
b = 3  # semi-minor axis
n = 2.5  # exponent (n > 2 for more rectangular, 1 < n < 2 for more oval-like)

# Create a meshgrid for the eta (angular coordinate), which ranges from 0 to 2pi
eta_vals = np.linspace(0, 2 * np.pi, 400)

# Calculate the corresponding xi values using the superellipse equation in elliptical coordinates
xi_vals = superellipse_elliptical(eta_vals, a, b, n)

# Plotting the result
plt.figure(figsize=(8, 8))
plt.plot(eta_vals, xi_vals, color='blue')
plt.title("Superellipse in Elliptical Coordinates")
plt.xlabel("Eta (Angular Coordinate)")
plt.ylabel("Xi (Radial Distance)")
plt.grid(True)
plt.show()
