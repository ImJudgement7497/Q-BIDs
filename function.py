import numpy as np

# def function(x, y, a=4):
#     theta = np.arctan2(y, x)  # Compute theta from Cartesian coordinates
#     r_curve = a * (3 * np.cos(theta) - np.cos(3 * theta))  # Radius from parametric equation
#     r_point = np.sqrt(x**2 + y**2)  # Radius of (x, y)
    
#     return r_point - r_curve  # Difference in radius (implicit function)

# def function(x, y):
#     
#     return (r-2)**2 + ((theta**2) / 4) - 1

r""" Nephroid Equation for constant a"""
# def function(x, y, a=1):
#     inside = x**2 + y**2 - 4*(a**2)
#     return (inside**3) - 108 * (y**2) * (a**4)

r""" Cardioid Equation for constant a"""
def function(x, y, a=2):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    
    return r - 2*a*(1-np.cos(theta))
    
    
# r"""Polar Nephroid Equation for constant a"""
# def function(x, y, a=1):
#     r = np.sqrt(x**2 + y**2)
#     theta = np.arctan2(y, x)
    
#     first = (r/(2*a))**(2/3)
#     second = np.abs(np.sin(theta/2)) ** (2/3)

#     third = np.abs(np.cos(theta/2)) ** (2/3)
    
#     return first - second - third

r"""Polar Nephroid Equation in terms of x and y"""
# def function(x, y, a=1):
#     first = np.abs(x/(2*a))**(2/3)
#     second = np.abs(np.sin(y/2)) ** (2/3)

#     third = np.abs(np.cos(y/2)) ** (2/3)
    
#     return first - second - third



# def function(x, y):
#     c = 10
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
    
#     return xi**2 + eta**2 - 1
