"""
	main.py

	Purpose: Script to run for visualisation.

	@author Thomas Caleb

	@version 2.0 12/12/2024
    
"""

import os, sys

from classes import Constants, SpacecraftParameters, Dataset
from plot_2d import plot_2d
from plot_2d_pdf import plot_2d_pdf
from plot_thrust_profile import plot_thrust_profile
from plot_nli import plot_nli
from plot_double_integrator import plot_double_integrator_u, plot_double_integrator_x
from plot_hypervolume import plot_hypervolume
from compute_LAM import compute_LAM
import matplotlib.pyplot as plt

global ALPHA_0_GMM, ALPHA_1_GMM
ALPHA_0_GMM = 0.5495506294920584 # Central weight [-]
ALPHA_1_GMM = 0.225224685253970 # Lateral weight [-]

"""
    Returns a dataset from a filename.

"""
def get_dataset(file_name):
    
    dataset = Dataset()
    dataset.read(file_name)
    return dataset


"""
    Function that runs first.

"""
if __name__ == "__main__":
    # Change directory to L-SODA for Windows
    # os.chdir("../../")
    
    # Get arguments.
    list_arguments = sys.argv[1:]
    test_case = int(list_arguments[0])
    
    LOADS_max_depth = float(list_arguments[1])
    inv_LOADS_max_depth = int(LOADS_max_depth)
    if LOADS_max_depth != 0:
        inv_LOADS_max_depth = int(1.0/float(list_arguments[1]))
    beta = float(list_arguments[2])
    inv_beta = int(beta)
    if beta != 0:
        inv_beta = int(1.0/float(list_arguments[2]))
    T2m_ratio = list_arguments[3]
    ToF = int(list_arguments[4])
    DDP_type = int(list_arguments[5])
    show_sample = int(list_arguments[6])
    
    extra = "" # At the end of the file name.
    
    # Make two names depending ontest case.
    file_name = "./data/"
    file_name_robust = file_name + "robust_trajectory/"
    file_name_sample = file_name + "sample_trajectory/"
    if test_case == 0:
        file_name_robust = file_name_robust + 'double_integrator'
        file_name_sample = file_name_sample + 'double_integrator'
    elif test_case == 1:
        file_name_robust = file_name_robust + 'tbp_SUN_lt_earth_to_mars'
        file_name_sample = file_name_sample + 'tbp_SUN_lt_earth_to_mars'
    elif test_case == 2:
        file_name_robust = file_name_robust + 'cr3bp_EARTH_MOON_lt_haloL2_to_haloL1'
        file_name_sample = file_name_sample + 'cr3bp_EARTH_MOON_lt_haloL2_to_haloL1'
    elif test_case == 3:
        file_name_robust = file_name_robust + 'cr3bp_EARTH_MOON_lt_nrho_to_dro'
        file_name_sample = file_name_sample + 'cr3bp_EARTH_MOON_lt_nrho_to_dro'
    elif test_case == 4:
        file_name_robust = file_name_robust + 'cr3bp_EARTH_MOON_lt_dro_to_dro' 
        file_name_sample = file_name_sample + 'cr3bp_EARTH_MOON_lt_dro_to_dro' 
    elif test_case == 5:
        file_name_robust = file_name_robust + 'cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2' 
        file_name_sample = file_name_sample + 'cr3bp_EARTH_MOON_lt_lyapunovL1_to_lyapunovL2' 
    
    if extra != "":
        extra = "_" + extra
    file_name_robust = file_name_robust + "_" + str(inv_LOADS_max_depth) + "_" + str(inv_beta) + "_" + T2m_ratio + "_" + str(ToF) + "_" + str(DDP_type) + extra + ".dat"
    file_name_sample = file_name_sample + "_" + str(inv_LOADS_max_depth) + "_" + str(inv_beta) + "_" + T2m_ratio + "_" + str(ToF) + "_" + str(DDP_type) + extra + ".dat"
    
    # Keywords to generate 2d plots and thurst plot.
    list_2d = ["tbp", "lyapunov", "dro", "leo", "meo", "halo"]
    
    # Load datasets.
    dataset_robust = get_dataset(file_name_robust) # Robust trajectory.
    dataset_robust.file_name = file_name_robust
    if show_sample == 1:
        dataset_sample = get_dataset(file_name_sample)
        dataset_sample.file_name = file_name_sample
    else:
        dataset_sample = Dataset()
    
    # Plots
    
    # Only for double integrator (test case 0).
    if ("double_integrator" in file_name_robust):
        plot_double_integrator_u(dataset_robust, dataset_sample)
        plot_double_integrator_x(dataset_robust, dataset_sample)

    # Astrodynamics test cases.
    else:   
        plot_thrust_profile(dataset_robust, dataset_sample)
        # compute_LAM(dataset_robust, dataset_sample)
        for i in list_2d:
            if i in file_name_robust:
                plot_2d(dataset_robust, dataset_sample)
                plot_2d_pdf(dataset_robust, dataset_sample)
                break
    