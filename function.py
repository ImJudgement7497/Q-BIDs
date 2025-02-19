import numpy as np

# def function(x, y):
#     r = np.sqrt(x**2 + y**2)
#     theta = np.arctan2(y, x)
#     return (r - 2)**2 + (theta**2) / 4 - 1

def function(x, y):
    c = 10
    B = x**2 + y**2 - c**2
    p = (-B + np.sqrt((B**2) + 4 * (c**2) * (y**2))) / (2 * (c**2))
    q = (-B - np.sqrt((B**2) + 4 * (c**2) * (y**2))) / (2 * (c**2))
    
    eta_0 = np.arcsin(np.sqrt(p))
    
    # Initialize eta with zeros
    eta = np.zeros_like(x)
    
    # Quadrant I (x >= 0, y >= 0)
    eta[(x >= 0) & (y >= 0)] = eta_0[(x >= 0) & (y >= 0)]
    
    # Quadrant II (x < 0, y >= 0)
    eta[(x < 0) & (y >= 0)] = np.pi - eta_0[(x < 0) & (y >= 0)]
    
    # Quadrant III (x =< 0, y < 0)
    eta[(x <= 0) & (y < 0)] = np.pi + eta_0[(x <= 0) & (y < 0)]
    
    # Quadrant IV (x > 0, y < 0)
    eta[(x > 0) & (y < 0)] = 2 * np.pi - eta_0[(x > 0) & (y < 0)]
    
    xi = 0.5 * np.log(1 - 2*q + 2 * np.sqrt((q**2) - q))
    
    return xi**2 + eta**2 - 1
