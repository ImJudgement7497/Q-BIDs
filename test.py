import numpy as np
import matplotlib.pyplot as plt

phi = np.linspace(0, 2*np.pi, 300)

a = 1
n = 2

quad = n**2 + 2*n + 2

r = np.sqrt(((a**2)/(n**2)) * (quad - 2 * (n+1) * np.cos(n * phi)))

top = (n+1) * np.sin(phi) - np.sin((n+1) * phi)
bot = (n+1) * np.cos(phi) - np.cos((n+1) * phi)
theta = np.arctan(top / bot)

# Convert to Cartesian coordinates for visualization in standard xy-plane
x = r * np.cos(theta)
y = r * np.sin(theta)

# Plot in polar coordinates
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(theta, r, label=r"$r = 2 + \cos(3t)$", color='b')

# Decorations
ax.set_title("Parametric Polar Plot")
ax.legend()
plt.savefig("parametric_polar_plot.png")

# Alternative: Plot in Cartesian coordinates
plt.figure()
plt.plot(x, y, label="Cartesian Projection", color="r")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")
plt.title("Cartesian Projection of Polar Curve")
plt.axhline(0, color='black', linewidth=1)
plt.axvline(0, color='black', linewidth=1)
plt.legend()
plt.grid()
plt.savefig("cartesian.png")