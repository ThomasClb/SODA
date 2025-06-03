"""
	plot_2d.py

	Purpose: Implements the functions to produce 2D plots of tranfers.

	@author Thomas Caleb

	@version 2.0 12/12/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import scipy.interpolate as interpolate
from scipy.stats import chi2, multivariate_normal

from misc import get_Lagrange_point
from classes import Dataset

plt.rcParams.update({'font.size': 15})
ALPHA_0_GMM = 0.5495506294920584 # Central weight [-]
ALPHA_1_GMM = 0.225224685253970 # Lateral weight [-]


"""
    Plots the specific points for a given dynamical system (two-body and 
    circular-restricted three-body problem).

"""
def plot_system_points(dataset, axis_0, axis_1, ax,
                       list_colors, list_markers,
                       list_plots, denormalise):    
    # Discriminate dynamical system
    if dataset.dynamical_system.startswith("TBP"):
        # Make points
        coord_center = np.array([0, 0, 0])
        
        # Denormalise
        if denormalise:
            LU = dataset.spacecraft_parameters.constants.LU
            coord_center = coord_center*LU
        
        # Plot
        if list_plots[0]: # Central body
            ax.scatter(coord_center[axis_0], coord_center[axis_1],
                       color=list_colors[0], marker=list_markers[0])
            
            # Print name
            central_body_name = dataset.dynamical_system.split(" ")[1]
            ax.text(coord_center[axis_0], coord_center[axis_1],
                    " " + central_body_name)

    elif dataset.dynamical_system.startswith("CR3BP"):
        # Make points
        coord_P1 = np.array([-dataset.spacecraft_parameters.constants.MU, 0, 0])
        coord_P2 = np.array([1-dataset.spacecraft_parameters.constants.MU, 0, 0])
        coord_L1 = get_Lagrange_point(dataset, 1)
        coord_L2 = get_Lagrange_point(dataset, 2)
        
        # Denormalise
        if denormalise:
            LU = dataset.spacecraft_parameters.constants.LU
            coord_P1 *= LU
            coord_P2 *= LU
            coord_L1 *= LU
            coord_L2 *= LU
            
        # Get names
        primary_names = dataset.dynamical_system.split(" ")[1].split("-")

        offset_1_x, offset_1_y, offset_2_x, offset_2_y, offset_p2_x, offset_p2_y = [0, 0, 0, 0, 0, 0]
        if ("halo" in dataset.file_name and axis_1 == 2):
            offset_1_x = 0.0
            offset_1_y = -0.025
            offset_2_x = -0.02
            offset_2_y = -0.025
        if ("halo" in dataset.file_name and axis_1 == 1):
            offset_1_x = -0.02
            offset_1_y = 0.025
            offset_2_x = -0.03
            offset_2_y = 0.025
        if ("nrho" in dataset.file_name):
            offset_1_x = 0
            offset_1_y = 0
            offset_2_x = -0.05
            offset_2_y = 0
            offset_p2_x = -0.11
        if ("dro_to_dro" in dataset.file_name):
            offset_1_x = 0
            offset_1_y = 0.02
            offset_2_x = -0.075
            offset_2_y = 0.02
            offset_p2_x = -0.1
            offset_p2_y = 0.02

        # Plots
        if list_plots[0]: # First primary
            ax.scatter(coord_P1[axis_0], coord_P1[axis_1],
                       color=list_colors[0], marker=list_markers[0])
            ax.text(coord_P1[axis_0], coord_P1[axis_1],
                    " " + primary_names[0])
            
        
        if list_plots[1]: # Second primary
            ax.scatter(coord_P2[axis_0], coord_P2[axis_1],
                       color=list_colors[1], marker=list_markers[1])
            ax.text(coord_P2[axis_0]+offset_p2_x, coord_P2[axis_1]+offset_p2_y,
                    " " + primary_names[1])
            
        if list_plots[2]: # Lagrange point 1
            ax.scatter(coord_L1[axis_0], coord_L1[axis_1],
                       color=list_colors[2], marker=list_markers[2])
            ax.text(coord_L1[axis_0]+offset_1_x, coord_L1[axis_1]+offset_1_y,
                    " $L_1$")
        
        if list_plots[3]: # Lagrange point 2
            ax.scatter(coord_L2[axis_0], coord_L2[axis_1],
                       color=list_colors[3], marker=list_markers[3])  
            ax.text(coord_L2[axis_0]+offset_2_x, coord_L2[axis_1]+offset_2_y,
                    " $L_2$")
        
"""
    Plots the departure and arrival points of a transfer.

"""    
def plot_departure_arrival(dataset, axis_0, axis_1, ax,
                           list_colors, list_markers,
                           denormalise):    
    # Retrieve data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i].copy()
        if (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i].copy()
                
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        data_state *= LU
            
    # Plot departure and arrival
    ax.scatter(data_departure[axis_0, 0], data_departure[axis_1, 0],
        s=10,
        color=list_colors[0], marker=list_markers[0])
    ax.text(data_departure[axis_0, 0], data_departure[axis_1, 0],
            " $x_0$",
              zorder=100)
    ax.scatter(data_arrival[axis_0, 0], data_arrival[axis_1, 0],
        s=10,
        color=list_colors[1], marker=list_markers[1])
    ax.text(data_arrival[axis_0, 0], data_arrival[axis_1, 0],
            " $x_t$",
              zorder=100)
    
"""
    Plots the thrust vectors along the trajectory.

"""
def plot_thrust_vector(dataset, axis_0, axis_1, ax,
    alpha, color, denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal control"):
            data_control = dataset.list_datasets[i].copy()
    coord_0 = data_state[axis_0,:]
    coord_1 = data_state[axis_1,:]
    ucoord_0 = data_control[axis_0,:]
    ucoord_1 = data_control[axis_1,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        coord_0 *= LU
        coord_1 *= LU
        ucoord_0 *= THRUSTU
        ucoord_1 *= THRUSTU
    
    # Plot arrows
    ax.quiver(coord_0[:-1], coord_1[:-1],
              ucoord_0, ucoord_1,
              color=color, label='Thrust',
              alpha=alpha,
              zorder=10)
    
"""
    Plots the departure and arrival orbits of a transfer.

"""
def plot_reference_orbits(dataset, axis_0, axis_1, ax,
                          alpha,
                          list_colors, list_linestyles,
                          denormalise):    
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i].copy()
    x_dep_0 = data_departure[axis_0,:]
    y_dep_1 = data_departure[axis_1,:]
    x_arr_0 = data_arrival[axis_0,:]
    y_arr_1 = data_arrival[axis_1,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        x_dep_0 *= LU
        y_dep_1 *= LU
        x_arr_0 *= LU
        y_arr_1 *= LU
    
    # Plot orbits
    ax.plot(x_dep_0, y_dep_1,
            label='Departure orbit',
            alpha=alpha,
            color=list_colors[0], linestyle=list_linestyles[0])
    ax.plot(x_arr_0, y_arr_1, 
            label='Arrival orbit',
            alpha=alpha,
            color=list_colors[1], linestyle=list_linestyles[1])


def plot_sample(dataset, dataset_sample, axis_0, axis_1, ax, plot_CL, max_sample_size,
                transparancy, color, linewidth,
                denormalise, interpolation, interpolation_rate):

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_coord_0 = []
    list_coord_1 = []
    list_alpha = []
    list_data_Sigma = []
    list_center = []
    list_inv_cov = []
    quad = np.zeros((2,2))
    for k in range(nb_GMM):
        for i in range(nb_datasets):
            if (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
                list_names_state = dataset.list_dataset_names[i].copy()
                list_coord_0.append(data_state[axis_0,:])
                list_coord_1.append(data_state[axis_1,:])
                list_center.append(np.array([data_state[axis_0,0], data_state[axis_1,0]]))
            elif (dataset.list_dataset_names[i][0] == "Sigma GMM " + str(k)):
                data_Sigma = dataset.list_datasets[i].copy()
                d = int(np.sqrt(len(data_Sigma[:,0])))
                N = len(data_Sigma[0,:])
                list_data_Sigma.append(data_Sigma)
                Sigma = data_Sigma[:,0]
                quad[0,0] = Sigma[axis_0 + d*axis_0]
                quad[1,1] = Sigma[axis_1 + d*axis_1]
                quad[0,1] = Sigma[axis_0 + d*axis_1]
                quad[1, 0] = quad[0,1]
                list_inv_cov.append(np.linalg.inv(quad))
            elif (dataset.list_dataset_names[i][0] == "Splitting history GMM " + str(k)):
                data_history = dataset.list_datasets[i].copy()
                alpha = 1
                if data_history.shape[1] != 0:
                    for j in range(len(data_history[0,:])):
                        if data_history[1,j] == 0:
                            alpha = alpha*ALPHA_0_GMM
                        elif abs(data_history[1,j]) == 1:
                            alpha = alpha*ALPHA_1_GMM
                list_alpha.append(alpha)

    # Compute mean
    coord_0 = list_alpha[0]*list_coord_0[0]
    coord_1 = list_alpha[0]*list_coord_1[0]
    for k in range(1, nb_GMM):
        coord_0 = coord_0 + list_alpha[k]*list_coord_0[k]
        coord_1 = coord_1 + list_alpha[k]*list_coord_1[k]

    # Get sample
    sample_size = min(max_sample_size, int(len(dataset_sample.list_dataset_names)/4))
    coord_0_sample = []
    coord_1_sample = []
    if sample_size != 0:
        nb_datasets = len(dataset_sample.list_dataset_names)
        for i in range(0, sample_size):
            for j in range(nb_datasets):
                if (dataset_sample.list_dataset_names[j][0] == "Sample state " + str(i)):
                    data_control_sample = dataset_sample.list_datasets[j].copy()
            coord_0_sample.append(data_control_sample[axis_0,:])
            coord_1_sample.append(data_control_sample[axis_1,:])

    # Normalisation
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        for k in range(nb_GMM):
            list_data_state[k][axis_0,:] *= LU
            list_data_state[k][axis_1,:] *= LU

        # Sample
        if sample_size != 0:
            for i in range(sample_size):
                coord_0_sample[i] = LU*coord_0_sample[i]
                coord_1_sample[i] = LU*coord_1_sample[i]
        
        # Labels
        list_names_state[axis_0 + 1] = list_names_state[axis_0 + 1].replace(
            "LU", "km")
        list_names_state[axis_1 + 1] = list_names_state[axis_1 + 1].replace(
            "LU", "km")

    if interpolation:
        t_old = np.linspace(0, 1, len(coord_0))  
        t_new = np.linspace(0, 1, interpolation_rate*len(coord_0))  

        # Offset
        if sample_size != 0:
            for i in range(sample_size):
                # Unpack
                coord_0_i = coord_0_sample[i]
                coord_1_i = coord_1_sample[i]
                coord_i = np.array([coord_0_i[0], coord_1_i[0]])

                # Find split
                list_maha = []
                for k in range(nb_GMM):
                    list_maha.append((coord_i - list_center[k]).T@list_inv_cov[k]@(coord_i - list_center[k]))
                index = np.argmin(list_maha)
                coord_0_split = list_coord_0[index]
                coord_1_split = list_coord_1[index]

                # Scale in split
                coord_0_i = coord_0_split + (coord_0_i - coord_0_split)
                coord_1_i = coord_1_split + (coord_1_i - coord_1_split)

                # Scale
                coord_0_sample[i] = coord_0 + (coord_0_i - coord_0)
                coord_1_sample[i] = coord_1 + (coord_1_i - coord_1)
        
        # Interpolate
        if sample_size != 0:
            for i in range(sample_size):
                coord_0_sample[i] = interpolate.interp1d(
                    t_old, coord_0_sample[i], kind='cubic')(t_new)
                coord_1_sample[i] = interpolate.interp1d(
                    t_old, coord_1_sample[i], kind='cubic')(t_new)

    # Sample
    if sample_size != 0:
        for i in range(sample_size):
            ax.plot(coord_0_sample[i], coord_1_sample[i],
                alpha=alpha,
                linewidth=linewidth,
                color=color,
                zorder=10)
    
"""
    Plots the thrust vectors along the trajectory.

"""
def plot_thrust_vector(dataset, axis_0, axis_1, ax,
    alpha, color, denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Nominal state GMM 0"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal control GMM 0"):
            data_control = dataset.list_datasets[i].copy()
    coord_0 = data_state[axis_0,:]
    coord_1 = data_state[axis_1,:]
    ucoord_0 = data_control[axis_0,:]
    ucoord_1 = data_control[axis_1,:]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        coord_0 *= LU
        coord_1 *= LU
        ucoord_0 *= THRUSTU
        ucoord_1 *= THRUSTU
    
    # Plot arrows
    ax.quiver(coord_0[:-1], coord_1[:-1],
              ucoord_0, ucoord_1,
              color=color, label='Thrust',
              alpha=alpha,
              zorder=10)

"""
    Plots a 2D plot from a given dataset.

"""
def plot_2d(dataset, dataset_sample=Dataset()):
    
    # Settings
    
    dpi = 200
    
    # Axes
    list_axis = [[0, 1]]
    if "halo" in dataset.file_name or "nrho" in dataset.file_name:
        list_axis = [[0, 1], [0, 2]]
    
    # Ellipses
    ellipse_alpha = 0.5
    ellipse_sampling = 7
    ellipse_nb_points = 1001
    ellipse_levels = np.array([1e-240, 1e-180, 1e-120, 1e-60, 1])
    plot_CL = True

    # System points
    if dataset.dynamical_system.startswith("CR3BP"):
        list_colors_system_points = ["black", "black", "black", "black"]
        list_markers_system_points = ["o", "s", "o", "o"]
        list_plots_system_points = [False, True, True, True]
        
    elif dataset.dynamical_system.startswith("TBP"):
        list_colors_system_points = ["black"]
        list_markers_system_points = ["o"]
        list_plots_system_points = [True]
    
    # Departure arrival
    list_colors_departure_arrival = ["black", "black"]
    list_markers_departure_arrival = ["^", "v"]
    
    # Reference orbits
    alpha_reference = 0.5
    list_colors_reference = ["#7f7f7f", "#7f7f7f"]
    list_linestyles_references = ["dotted", "dashed"]

    # Sample
    color_sample = "#9a7bb5"
    max_sample_size = 200
    sample_alpha = 0.3
    sample_linewidth = 0.5

    # Trajectory
    color_trajectory = "black"
    linewidth_trajectory = 2
    
    # Thrust
    thrust_color = "#e63946"
    thurst_alpha = 0.8
    show_thrust = True

    # Normalisation
    denormalise = False

    # Interpolation
    interpolation = True
    interpolation_rate = 15

    # Legend
    show_legend = False
    legend_loc = "lower right"
    
    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)
    max_trust = dataset.spacecraft_parameters.thrust

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_data_state = []
    list_alpha = []
    for k in range(nb_GMM):
        for i in range(nb_datasets):
            if (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
                list_names_state = dataset.list_dataset_names[i].copy()
                list_data_state.append(data_state)
                if k==0:
                    date_state_0 = data_state

    for axies in list_axis:
        axis_0 = axies[0]
        axis_1 = axies[1]

        coord_0 = date_state_0[axis_0,:]
        coord_1 = date_state_0[axis_1,:]

        # Normalisation
        if denormalise:        
            # Labels
            list_names_state[axis_0 + 1] = list_names_state[axis_0 + 1].replace(
                "LU", "km")
            list_names_state[axis_1 + 1] = list_names_state[axis_1 + 1].replace(
                "LU", "km")

        # Create plot
        fig = plt.figure(dpi=dpi)
        ax = fig.add_subplot()
        ax.set_aspect("equal")

        # Set labels
        ax.set_xlabel(list_names_state[axis_0 + 1])
        ax.set_ylabel(list_names_state[axis_1 + 1])

        plot_sample(dataset, dataset_sample, axis_0, axis_1, ax, plot_CL, max_sample_size,
                    sample_alpha, color_sample, sample_linewidth,
                    denormalise, interpolation, interpolation_rate)
        
        # Plot reference orbits
        plot_reference_orbits(dataset, axis_0, axis_1, ax,
                              alpha_reference,
                              list_colors_reference,
                              list_linestyles_references,
                              denormalise)
        
        # Plot departure and arrival points
        plot_departure_arrival(dataset, axis_0, axis_1, ax,
                               list_colors_departure_arrival,
                               list_markers_departure_arrival,
                               denormalise)

        # Plot Thrust
        if show_thrust:
            plot_thrust_vector(dataset, axis_0, axis_1,
                               ax, thurst_alpha, thrust_color,
                               denormalise)
        
        # Plot system points
        plot_system_points(dataset, axis_0, axis_1, ax,
                           list_colors_system_points,
                           list_markers_system_points,
                           list_plots_system_points,
                           denormalise)

        # Plot trajectory
        ax.plot(coord_0, coord_1,
                color=color_trajectory,
                label='Trajectory',
                zorder=11,
                linewidth=linewidth_trajectory)

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
            signature = ("_2d_" + str(axis_0) + "_" + str(axis_1)
                + "." + saving_format)
            file_name = file_name.replace(
                ".dat", signature)
            
            # Save
            plt.savefig(file_name, bbox_inches='tight')    
           
        if show_plot:   
            plt.show()
        else:
            plt.close(fig)
    