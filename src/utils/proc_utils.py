import numpy as np

# *******************************************************************

def gaussian_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    '''
    This function is used to fit a 2D Gaussian to the solar scan image.
    This function is now (2025-03-31) obsolete. We now use the 2D paraboloid.
    '''
    x, y = coords
    return amplitude * np.exp(-(((x-x0)**2)/(2*sigma_x**2) + ((y-y0)**2)/(2*sigma_y**2)))

#**************************************************************************

def paraboloid_2d(coords, amplitude, x0, y0, sigma_x, sigma_y):
    '''
    This function is used to fit a 2D paraboloid to the solar scan image.
    '''
    x, y = coords
    return amplitude - 4*np.log(2) * ((x-x0)**2/sigma_x**2 + (y-y0)**2/sigma_y**2)

#**************************************************************************