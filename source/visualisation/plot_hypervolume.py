"""
	plot_cdf.py

	Purpose: Implements the functions to plot the ellipsoïd volume.

	@author Thomas Caleb

	@version 1.0 12/12/2024
    
"""

import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy.stats import chi2, multivariate_normal

from classes import Dataset

# Plots the ellipsoïd volume for a given trajectory.
def plot_hypervolume(dataset):
    # Settings
    
    dpi = 200

    # Main data
    state_color = "black"
    state_linewidth = 2   
    control_color = "#e63946"
    control_linewidth = 2

    # Normalisation
    denormalise = True
    
    # Legend
    show_legend = True
    legend_loc = "lower left"

    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False

     # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Sigma"):
            data_Sigma = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal control"):
            data_control = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Feedback gain"):
            data_gains = dataset.list_datasets[i].copy()
    dt = data_state[7,0]
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    u = np.sqrt(ux*ux + uy*uy + uz*uz)
    N = len(data_state[7,:])
    x_label = "Time [TU]"
    y_label = "Uncertainty hypervolume [-]"

    # Get beta
    list_file_name = dataset.file_name.split("_")
    inv_beta = float(list_file_name[-4])
    beta = 0.05
    if inv_beta != 0:
        beta = 1/inv_beta

    # Normalisation
    if denormalise:
        TU = dataset.spacecraft_parameters.constants.TU
        dt *= (TU/86400)
        
        # Labels
        x_label = x_label.replace("TU", "days")
        
    # Compute confidence intervals.
    
    # Get state hypervoume.
    list_hypervolume_state = []
    for i in range(N):
        if i==0:
            hypervolume_1 = np.linalg.det(np.linalg.cholesky(data_Sigma[:,0].reshape((8,8))[:6, :6]))

        # From https://fr.wikipedia.org/wiki/Ellipso%C3%AFde#Volume
        Sigma_x = data_Sigma[:,i].reshape((8,8))
        hyperellipse_volume = (np.linalg.det(np.linalg.cholesky(Sigma_x[:6, :6]))/hypervolume_1)**(1/6.0)
        list_hypervolume_state.append(hyperellipse_volume)

    # Get control hypervolume.
    list_hypervolume_control = []
    for i in range(N-1):
       
        # Build Sigma.
        K = data_gains[:,i].reshape((3,8))
        Sigma_x = data_Sigma[:,i].reshape((8,8))
        Sigma = K@Sigma_x@K.T

        # From https://fr.wikipedia.org/wiki/Ellipso%C3%AFde#Volume
        if i==0:
            hypervolume_1 = np.linalg.det(np.linalg.cholesky(Sigma))
        hyperellipse_volume = (np.linalg.det(np.linalg.cholesky(Sigma))/hypervolume_1)**(1/3.0) # ADD
        list_hypervolume_control.append(hyperellipse_volume)
    
    # Create plot
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot()

    # Set labelss
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_yscale('log')

    # Plots
    t = dt*np.array(range(N))
    ax.plot(dt*np.array(range(N)), list_hypervolume_state,
        label="State", color=state_color, linewidth=state_linewidth)
    ax.plot(dt*np.array(range(N-1)), list_hypervolume_control,
        label="Control", color=control_color, linewidth=control_linewidth)

    # Layout
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        plt.grid(alpha=0.5, color="#d3d3d3")

    if show_legend: 
        plt.legend(loc=legend_loc, frameon=True, fancybox=True, facecolor='white', edgecolor='black')
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "robust_trajectory", "plots")
                
        # Add signature
        signature = ("_hypervolume" + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
