"""
    plot_thrust_profile.py

    Purpose: Implements the functions to produce thrust profiles of tranfers.

    @author Thomas Caleb

    @version 2.0 12/12/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2, multivariate_normal

from classes import Dataset

plt.rcParams.update({'font.size': 15})

"""
    Turns a continuous vector into a stairs shape larger vector.

"""
def make_stairs(dt, u):
    nb_point_stair = 50
    nb_point = len(dt)
    t = []
    u_stairs = []
    current_time = 0
    for i in range(nb_point):
        dt_i = dt[i]/(float(nb_point_stair))
        u_i = u[i]
        for j in range(nb_point_stair):
            t.append(current_time)
            u_stairs.append(u_i)
            current_time = current_time + dt_i
    return [np.array(t), np.array(u_stairs)]
    

"""
    Plots a thrust profile with a stairs shape.

"""
def plot_stairs_profile(dt, u, ax, label, color,
                        linestyle="solid"):
    t, u_stairs = make_stairs(dt, u) 
    ax.plot(t, u_stairs, label=label,
            color=color, linestyle=linestyle)

"""
    Plots a the expected confidence interval for a profile.

"""
def plot_control_distribution(dataset, ax, scale, color_line, denormalise):
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
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    dt = data_state[7,:-1]
    u = np.sqrt(ux*ux + uy*uy + uz*uz)
    u_down = 0 + u
    u_up = 0 + u
    max_trust = dataset.spacecraft_parameters.thrust
    x_label = "Time [TU]"
    y_label = "Thrust norm [THRUSTU]"

    # Get beta
    list_file_name = dataset.file_name.split("_")
    inv_beta = float(list_file_name[-4])
    beta = 0.05
    if inv_beta != 0:
        beta = 1/inv_beta

    # Normalisation
    if denormalise:
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        TU = dataset.spacecraft_parameters.constants.TU
        ux *= THRUSTU
        uy *= THRUSTU
        uz *= THRUSTU
        u *= THRUSTU
        max_trust *= THRUSTU
        dt *= (TU/86400)
        
        # Labels
        x_label = x_label.replace("TU", "days")
        y_label = y_label.replace("THRUSTU", "N")
        
    # Compute confidence intervals
    
    # Get Sigma
    for i in range(len(u)):
        K = data_gains[:,i].reshape((3,8))
        Sigma_x = data_Sigma[:,i].reshape((8,8))
        A = np.array([ux[i], uy[i], uz[i]]).T@K
        if u[i] != 0:
            A = A/u[i]
        Sigma = A@Sigma_x@A.T
        
        # Normalisation
        if denormalise:
            Sigma *= THRUSTU*THRUSTU
        
        # Scale
        bound = scale*np.sqrt(Sigma*chi2.ppf(1-beta, 2)) # First order transcription.
        u_up[i] = u[i] + bound
        u_down[i] = u[i] - bound    
        
    # Plot CI.
    t, u_stairs_up = make_stairs(dt, u_up)
    t, u_stairs_down = make_stairs(dt, u_down)
    ax.plot(t, u_stairs_up, label="Test",
            color=color_line, linestyle="dashed")
    ax.plot(t, u_stairs_down, label="Test",
            color=color_line, linestyle="dashed")

"""
    Plots a thrust pofile for a given transfer dataset.

"""
def plot_thrust_profile(dataset, dataset_sample=Dataset()):
    # Settings
    
    dpi = 200
    
    # Thrust norm
    color_thrust_norm = "black"
    label_thrust_norm  = "Thrust"
    linewidth_thrust_norm = 1.5
     
    # Max thrust norm
    color_thrust_max = "#e63946"
    linestyle_thrust_max = "dotted"
    label_thrust_max  = "Max thrust"
    
    # Robust thrust norm
    color_thrust_robust = "#4a90e2"
    linestyle_thrust_robust = "dashed"
    alpha_thrust_robust = 0.8
    linewidth_thrust_robust = 0.8
    label_thrust_robust  = "Margins"
    if "mars" in dataset.file_name:
        scale_thrust_robust = 200
    elif "halo" in dataset.file_name:
        scale_thrust_robust = 2000
    elif "nrho" in dataset.file_name:
        scale_thrust_robust = 1000
    elif "dro_to_dro" in dataset.file_name:
        scale_thrust_robust = 200
    elif "lyapunov" in dataset.file_name:
        scale_thrust_robust = 3000

    # Sample norm
    sample_color = "#9a7bb5"
    sample_alpha = 0.3
    sample_linewidth = 0.5
    max_sample_size = 100
    
    # Normalisation
    denormalise = True
    
    # Legend
    show_legend = False
    legend_loc = "lower right"

    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False
    
    # Retreive data
    nb_datadets = len(dataset.list_dataset_names)
    for i in range(nb_datadets):
        if (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal control"):
            data_control = dataset.list_datasets[i].copy()
    ux = data_control[0,:]
    uy = data_control[1,:]
    uz = data_control[2,:]
    dt = data_state[7,:-1]
    u = np.sqrt(ux*ux + uy*uy + uz*uz)
    max_trust = dataset.spacecraft_parameters.thrust
    x_label = "Time [TU]"
    y_label = "Thrust norm [THRUSTU]"

    # Get sample
    sample_size = min(max_sample_size, int(len(dataset_sample.list_dataset_names)/4))
    ux_sample = []
    uy_sample = []
    uz_sample = []
    u_sample = []
    if sample_size != 0:
        nb_datadets = len(dataset_sample.list_dataset_names)
        for i in range(0, sample_size):
            for j in range(nb_datadets):
                if (dataset_sample.list_dataset_names[j][0] == "Sample control " + str(i)):
                    data_control_sample = dataset_sample.list_datasets[j].copy()
            ux_sample.append(data_control_sample[0,:])
            uy_sample.append(data_control_sample[1,:])
            uz_sample.append(data_control_sample[2,:])
            u_sample.append(
                np.sqrt(data_control_sample[0,:]*data_control_sample[0,:]
                    + data_control_sample[1,:]*data_control_sample[1,:]
                    + data_control_sample[2,:]*data_control_sample[2,:]))
    
    # Normalisation
    if denormalise:
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        TU = dataset.spacecraft_parameters.constants.TU
        ux *= THRUSTU
        uy *= THRUSTU
        uz *= THRUSTU
        u *= THRUSTU
        max_trust *= THRUSTU
        dt *= (TU/86400)

        # Sample
        if sample_size != 0:
            for i in range(sample_size):
                u_sample[i] = THRUSTU*u_sample[i]
                ux_sample[i] = THRUSTU*ux_sample[i]
                uy_sample[i] = THRUSTU*uy_sample[i]
                uz_sample[i] = THRUSTU*uz_sample[i]
        
        # Labels
        x_label = x_label.replace("TU", "days")
        y_label = y_label.replace("THRUSTU", "N")


    # Create plot
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot()

    # Set labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    # Axies limites
    plt.ylim((0 - 0.05*max_trust, 1.05*max_trust))
    
    # Plot sample
    if sample_size != 0:
        for i in range(sample_size):
            t, u_stairs = make_stairs(dt, u + scale_thrust_robust*(u_sample[i] - u))
            ax.plot(t, u_stairs, linewidth=sample_linewidth,
                alpha=sample_alpha,
                color=sample_color)

    # Plot errors.
    plot_control_distribution(dataset, ax,
                              scale_thrust_robust, color_thrust_robust, denormalise)

    # Plot Thrust
    t, u_stairs = make_stairs(dt, u)
    ax.plot(t, u_stairs, label=label_thrust_norm,
        linewidth=linewidth_thrust_norm,
        color=color_thrust_norm)
       
    # Plot max thrust
    t, u_stairs = make_stairs(dt, u*0 + max_trust)
    ax.plot(t, u_stairs, label=label_thrust_max,
            color=color_thrust_max, linestyle=linestyle_thrust_max)

    # Layout
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        plt.grid(alpha=0.5, color="#d3d3d3")
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "robust_trajectory", "plots")
                
        # Add signature
        signature = ("_robust_thrust" + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    