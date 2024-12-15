"""
	plot_nli.py

	Purpose: Implements the functions to plot transfer's non-linearity index (NLI).

	@author Thomas Caleb

	@version 1.0 21/11/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt

from classes import Dataset
    

"""
    Plots the NLI for a given transfer dataset.
    https://doi.org/10.2514/1.G007271

"""
def plot_nli(dataset):
    # Settings
    
    dpi = 200
    
    # Thrust norm
    color_nli = "black"
     
    # Output
    denormalise = True
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False
    
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "NLI"):
            nli = dataset.list_datasets[i].copy()
        if (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
    dt = data_state[7,:]
    x_label = "Time [TU]"
    y_label = "NLI [-]"
    
    # Normalisation
    if denormalise:
        TU = dataset.spacecraft_parameters.constants.TU
        dt *= (TU/86400)
        
        # Labels
        x_label = x_label.replace("TU", "days")

    # Create plot
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot()

    # Set labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_yscale('log')
    
    # Make t
    t = np.zeros(nli.shape) + dt[0]
    for i in range(nli.shape[1]):
        if i > 0:
            t[:,i] = t[:,i - 1] + dt[i]

    ax.plot(t[0,:], nli[0,:], color=color_nli)

    # Layout
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        plt.grid(alpha=0.5, color="#d3d3d3")

    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "robust_trajectory", "plots")
                
        # Add signature
        signature = ("_robust_nli" + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    