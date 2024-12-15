"""
	plot_double_integrator.py

	Purpose: Implements the functions to produce plot for double_integrator.

	@author Thomas Caleb

	@version 2.0 12/12/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2, multivariate_normal

from misc import get_Lagrange_point
from classes import Dataset

"""
    Plots the excpected control errors.

"""
def plot_control_errors(dataset,ax, scale, list_colors):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Sigma"):
            data_Sigma = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal control"):
            data_control = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Feedback gain"):
            data_gains = dataset.list_datasets[i].copy()
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    ux_up = ux + 0
    ux_down = ux + 0
    uy_up = uy + 0
    uy_down = uy + 0   
    uz_up = uz + 0
    uz_down = uz + 0
    N = len(ux)
    t = range(N)

    # Get beta
    list_file_name = dataset.file_name.split("_")
    inv_beta = float(list_file_name[-4])
    beta = 0.05
    if inv_beta != 0:
        beta = 1/inv_beta

     # Get Sigma
    for i in range(len(ux)):
        K = data_gains[:,i].reshape((3,7))
        Sigma_x = data_Sigma[:,i].reshape((7,7))
        Sigma_k = K@Sigma_x@K.T
        
        # Scale
        bound = np.sqrt(chi2.ppf(1-beta, 3))

        # Compute bounds
        ux_up[i] += scale*bound*np.sqrt(Sigma_k[0, 0])
        ux_down[i] -= scale*bound*np.sqrt(Sigma_k[0, 0])
        uy_up[i] += scale*bound*np.sqrt(Sigma_k[1, 1])
        uy_down[i] -= scale*bound*np.sqrt(Sigma_k[1, 1])
        uz_up[i] += scale*bound*np.sqrt(Sigma_k[2, 2])
        uz_down[i] -= scale*bound*np.sqrt(Sigma_k[2, 2])

    # Plots
    ax.plot(t, ux_up,
        color=list_colors[0], linestyle="dashed", linewidth=0.5)
    ax.plot(t, ux_down,
        color=list_colors[0], linestyle="dashed", linewidth=0.5)
    ax.plot(t, uy_up,
        color=list_colors[1], linestyle="dashed", linewidth=0.5)
    ax.plot(t, uy_down,
        color=list_colors[1], linestyle="dashed", linewidth=0.5)
    ax.plot(t, uz_up,
        color=list_colors[2], linestyle="dashed", linewidth=0.5)
    ax.plot(t, uz_down,
        color=list_colors[2], linestyle="dashed", linewidth=0.5)

    return

"""
    Plots the excpected state errors.

"""
def plot_state_errors(dataset,ax, scale, list_colors):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Sigma"):
            data_Sigma = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
    x = data_state[0,:]
    y = data_state[1,:]
    z = data_state[2,:]
    x_up = x + 0
    x_down = x + 0
    y_up = y + 0
    y_down = y + 0   
    z_up = z + 0
    z_down = z + 0
    N = len(z)
    t = range(N)

    # Get beta
    list_file_name = dataset.file_name.split("_")
    beta = 1/float(list_file_name[-4])

     # Get Sigma
    for i in range(len(x)):
        Sigma_x = data_Sigma[:,i].reshape((7,7))
        
        # Scale
        bound = np.sqrt(chi2.ppf(1-beta, 6))

        # Compute bounds
        x_up[i] += scale*bound*np.sqrt(Sigma_x[0, 0])
        x_down[i] -= scale*bound*np.sqrt(Sigma_x[0, 0])
        y_up[i] += scale*bound*np.sqrt(Sigma_x[1, 1])
        y_down[i] -= scale*bound*np.sqrt(Sigma_x[1, 1])
        z_up[i] += scale*bound*np.sqrt(Sigma_x[2, 2])
        z_down[i] -= scale*bound*np.sqrt(Sigma_x[2, 2])

    # Plots.
    ax.plot(t, x_up,
        color=list_colors[0], linestyle="dashed", linewidth=0.5)
    ax.plot(t, x_down,
        color=list_colors[0], linestyle="dashed", linewidth=0.5)
    ax.plot(t, y_up,
        color=list_colors[1], linestyle="dashed", linewidth=0.5)
    ax.plot(t, y_down,
        color=list_colors[1], linestyle="dashed", linewidth=0.5)
    ax.plot(t, z_up,
        color=list_colors[2], linestyle="dashed", linewidth=0.5)
    ax.plot(t, z_down,
        color=list_colors[2], linestyle="dashed", linewidth=0.5)
    
    return

"""
    Plots the control for a given dataset.

"""
def plot_double_integrator_u(dataset, dataset_sample=Dataset()):
    
    # Settings
    dpi = 200
        
    # Departure arrival
    list_colors = ["blue", "green", "purple"]
    list_markers = ["o", "s", "^"]

    # Sample
    max_sample_size = 200
    sample_alpha = 0.3
    sample_linewidth = 0.1

    # Legend
    show_legend = True
    legend_loc_control = "upper left"
    legend_loc_state = "lower left"
    
    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Nominal control"):
            data_control = dataset.list_datasets[i].copy()
            list_names_control = dataset.list_dataset_names[i].copy()
    control_0 = data_control[0,:]
    control_1 = data_control[1,:]
    control_2 = data_control[2,:]
    N = len(control_1)

    # Get sample
    sample_size = min(max_sample_size, int(len(dataset_sample.list_dataset_names)/4))
    control_0_sample = []
    control_1_sample = []
    control_2_sample = []
    if sample_size != 0:
        nb_datadets = len(dataset_sample.list_dataset_names)
        for i in range(0, sample_size):
            for j in range(nb_datadets):
                if (dataset_sample.list_dataset_names[j][0] == "Sample control " + str(i)):
                    data_control_sample = dataset_sample.list_datasets[j].copy()
            control_0_sample.append(data_control_sample[0,:])
            control_1_sample.append(data_control_sample[1,:])
            control_2_sample.append(data_control_sample[2,:])


    # Change labels
    list_names_control[0 + 1] = list_names_control[0 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[1 + 1] = list_names_control[1 + 1].replace(
        "[THRUSTU]", "")
    list_names_control[2 + 1] = list_names_control[2 + 1].replace(
        "[THRUSTU]", "")

    # Create plot
    fig, ax1 = plt.subplots(1, 1, dpi=dpi)

    # Set labels
    ax1.set_xlabel("Stage [-]")
    ax1.set_ylabel("Controls [-]")
    
    # Plot control
    t = range(N)
    ax1.scatter(t, control_0,
            color=list_colors[0],
            marker=list_markers[0],
            label=list_names_control[1])
    ax1.scatter(t, control_1,
            color=list_colors[1],
            marker=list_markers[1],
            label=list_names_control[2][1:])
    ax1.scatter(t, control_2,
            color=list_colors[2],
            marker=list_markers[2],
            label=list_names_control[3][1:])
    plot_control_errors(dataset, ax1, 1, list_colors)

    # Sample
    if sample_size != 0:
        for i in range(sample_size):
            ax1.plot(t, control_0_sample[i],
                color=list_colors[0],
                alpha=sample_alpha,
                linewidth=sample_linewidth,
                zorder=1)
            ax1.plot(t, control_1_sample[i],
                color=list_colors[1],
                alpha=sample_alpha,
                linewidth=sample_linewidth,
                zorder=1)
            ax1.plot(t, control_2_sample[i],
                color=list_colors[2],
                alpha=sample_alpha,
                linewidth=sample_linewidth,
                zorder=1)

    # Layout
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        ax1.grid(alpha=0.5, color="#d3d3d3")
    
    if show_legend: 
        ax1.legend(loc=legend_loc_control)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "robust_trajectory", "plots")
                
        # Add signature
        signature = ("_u." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)

"""
    Plots states from a given dataset.

"""
def plot_double_integrator_x(dataset, dataset_sample=Dataset()):
    
    # Settings
    dpi = 200

    # Sample
    max_sample_size = 200
    sample_alpha = 0.3
    sample_linewidth = 0.1
        
    # Departure arrival
    list_colors = ["blue", "green", "purple"]
    list_markers = ["o", "s", "^"]

    # Legend
    show_legend = True
    legend_loc_control = "upper left"
    legend_loc_state = "lower left"
    
    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
            list_names_state = dataset.list_dataset_names[i].copy()
    coord_0 = data_state[0,:]
    coord_1 = data_state[1,:]
    coord_2 = data_state[2,:]
    N = len(coord_0) - 1

    # Get sample
    sample_size = min(max_sample_size, int(len(dataset_sample.list_dataset_names)/4))
    coord_0_sample = []
    coord_1_sample = []
    coord_2_sample = []
    if sample_size != 0:
        nb_datadets = len(dataset_sample.list_dataset_names)
        for i in range(0, sample_size):
            for j in range(nb_datadets):
                if (dataset_sample.list_dataset_names[j][0] == "Sample state " + str(i)):
                    data_control_sample = dataset_sample.list_datasets[j].copy()
            coord_0_sample.append(data_control_sample[0,:])
            coord_1_sample.append(data_control_sample[1,:])
            coord_2_sample.append(data_control_sample[2,:])

    # Change labels
    list_names_state[0 + 1] = list_names_state[0 + 1].replace(
        "[LU]", "")
    list_names_state[1 + 1] = list_names_state[1 + 1].replace(
        "[LU]", "")
    list_names_state[2 + 1] = list_names_state[2 + 1].replace(
        "[LU]", "")


    # Create plot
    fig, ax2 = plt.subplots(1, 1, dpi=dpi)

    # Set labels
    ax2.set_xlabel("Stage [-]")
    ax2.set_ylabel("States [-]")

    # Plot state
    t = range(N+1)
    ax2.scatter(t, coord_0,
            color=list_colors[0],
            marker=list_markers[0],
            label=list_names_state[1])
    ax2.scatter(t, coord_1,
            color=list_colors[1],
            marker=list_markers[1],
            label=list_names_state[2][1:])
    ax2.scatter(t, coord_2,
            color=list_colors[2],
            marker=list_markers[2],
            label=list_names_state[3][1:])
    plot_state_errors(dataset, ax2, 1, list_colors)

    # Sample
    if sample_size != 0:
        for i in range(sample_size):
            ax2.plot(t, coord_0_sample[i],
                color=list_colors[0],
                alpha=sample_alpha,
                linewidth=sample_linewidth,
                zorder=1)
            ax2.plot(t, coord_1_sample[i],
                color=list_colors[1],
                alpha=sample_alpha,
                linewidth=sample_linewidth,
                zorder=1)
            ax2.plot(t, coord_2_sample[i],
                color=list_colors[2],
                alpha=sample_alpha,
                linewidth=sample_linewidth,
                zorder=1)


    # Layout
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        ax2.grid(alpha=0.5, color="#d3d3d3")
    
    if show_legend: 
        ax2.legend(loc=legend_loc_state)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "robust_trajectory", "plots")
                
        # Add signature
        signature = ("_x." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    