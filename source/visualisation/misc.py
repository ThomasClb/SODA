"""
	misc.py

	Purpose: Implementation of varied functions needed for visualisation.

	@author Thomas Caleb

	@version 1.0 16/01/2024
    
"""

import numpy as np

from classes import Dataset


"""
    Computes the approximate position of the Lagrangian points
    Up to the order 4 in eps.
    https://fr.wikipedia.org/wiki/Point_de_Lagrange#Calcul_de_la_position_des_points_de_Lagrange

"""
def get_Lagrange_point(dataset, point):
    q = dataset.spacecraft_parameters.constants.MU
    eps = (q/3.0)**(1/3.0)
    
    # Discriminate for each Lagrange point
    if point == 1:  
        x_L = 1.0-q - eps + eps**2/3 + eps**3/9
    elif point == 2:  
        x_L = 1.0 - q + eps + eps**2/3 - eps**3/9
    elif point == 3:  
        x_L = -1.0 - 5*q/12 + + 1127*q**3/20736
    
    return np.array([x_L, 0.0, 0.0])
