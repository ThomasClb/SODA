"""
	plot_2d.py

	Purpose: Implements the functions to produce 2D plots of tranfers.

	@author Thomas Caleb

	@version 2.0 12/12/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from scipy.stats import chi2, multivariate_normal

from misc import get_Lagrange_point
from classes import Dataset

plt.rcParams.update({'font.size': 15})



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
        if (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
            
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        data_state *= LU
            
    # Plot departure and arrival
    ax.scatter(data_state[axis_0, 0], data_state[axis_1, 0],
               color=list_colors[0], marker=list_markers[0])
    ax.text(data_state[axis_0, 0], data_state[axis_1, 0],
            " $x_0$",
              zorder=100)
    ax.scatter(data_state[axis_0, -1], data_state[axis_1, -1],
               color=list_colors[1], marker=list_markers[1])
    ax.text(data_state[axis_0, -1], data_state[axis_1, -1],
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
            color=list_colors[0], linestyle=list_linestyles[0])
    ax.plot(x_arr_0, y_arr_1, 
            label='Arrival orbit',
            color=list_colors[1], linestyle=list_linestyles[1])


"""
    Plots the projections of the uncertainty ellipsoïds of a transfer.

"""
def plot_state_distribution(dataset, axis_0, axis_1, ax, plot_CL,
                            scale, alpha, color, linewidth, denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Sigma"):
            data_Sigma = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Der dynamics"):
            data_der = dataset.list_datasets[i].copy()
    LU = dataset.spacecraft_parameters.constants.LU

    # Get beta
    list_file_name = dataset.file_name.split("_")
    inv_beta = float(list_file_name[-4])
    beta = 0.05
    if inv_beta != 0:
        beta = 1/inv_beta
    
    # Plot ellispses
    quad = np.zeros((2,2))
    t = np.linspace(0, 2*np.pi, 50)
    d = int(np.sqrt(len(data_Sigma[:,0])))
    N = len(data_Sigma[0,:])
    bound = np.sqrt(chi2.ppf(1-beta, d - 2))
    if plot_CL: # Closed loop
        for i in range(N):
            # Get projected cov
            quad[0,0] = data_Sigma[axis_0 + d*axis_0, i]
            quad[1,1] = data_Sigma[axis_1 + d*axis_1, i]
            quad[0,1] = data_Sigma[axis_0 + d*axis_1, i]
            quad[1, 0] = quad[0,1]
            center = data_state[:,i]
            
            # Denormalise
            if denormalise:
                center *= LU
                quad *= LU*LU
            
            # Get eigenvalues
            eig, eigvectors = np.linalg.eig(quad)
           
            # Make ellipse shape
            Ell = bound*scale*np.array([np.sqrt(eig[0])*np.cos(t) , np.sqrt(eig[1])*np.sin(t)])  
           	
            # Rotation
            R_rot = np.array([eigvectors[0], eigvectors[1]])  
            Ell_rot = np.zeros((2,Ell.shape[1]))
            for i in range(Ell.shape[1]):
                Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
           
            # Plot
            plt.plot(center[axis_0]+Ell_rot[0,:], center[axis_1]+Ell_rot[1,:],
                     color=color, linewidth=linewidth, linestyle="dashed", alpha=alpha,
                     zorder=0)
    else: # Open loop
        for i in range(0,N):
            if i == 0:
                Sigma_i = data_Sigma[:, 0].reshape((d, d))

            else:
                der = data_der[:,i].reshape((d,d+3))[:d,:d]
                Sigma_i = der.T@Sigma_im1@der

            # Get projected cov
            quad[0,0] = Sigma_i[axis_0, axis_0]
            quad[1,1] = Sigma_i[axis_1, axis_1]
            quad[0,1] = Sigma_i[axis_0, axis_1]
            quad[1, 0] = quad[0,1]
            center = data_state[:,i]
            Sigma_im1 = Sigma_i
            
            # Denormalise
            if denormalise:
                center *= LU
                quad *= LU*LU
            
            # Get eigenvalues
            eig, eigvectors = np.linalg.eig(quad)
           
            # Make ellipse shape
            Ell = bound*scale*np.array([np.sqrt(eig[0])*np.cos(t) , np.sqrt(eig[1])*np.sin(t)])  
            
            # Rotation
            R_rot = np.array([eigvectors[0], eigvectors[1]])  
            Ell_rot = np.zeros((2,Ell.shape[1]))
            for i in range(Ell.shape[1]):
                Ell_rot[:,i] = np.dot(R_rot,Ell[:,i])
           
            # Plot
            plt.plot(center[axis_0]+Ell_rot[0,:], center[axis_1]+Ell_rot[1,:],
                     color=color_line, linewidth=1.0, linestyle="dashed", alpha=0.4,
                     zorder=0)
"""
    Plots a 2D plot from a given dataset.

"""
def plot_2d(dataset, dataset_sample=Dataset()):
    
    # Settings
    
    dpi = 200
    
    # Axes
    list_axis = [[0, 1]]
    if "halo" in dataset.file_name:
        list_axis = [[0, 1], [0, 2]]
    
    # Ellipses
    ellipse_color = "#4a90e2"
    ellipse_alpha = 0.7 # 0.3
    ellipse_linewidth = 0.8
    plot_CL = True
    if "mars" in dataset.file_name:
        ellipse_scale = 5e4
    elif "halo" in dataset.file_name:
        ellipse_scale = 3e5
    elif "nrho" in dataset.file_name:
        ellipse_scale = 2e5
    elif "dro_to_dro" in dataset.file_name:
        ellipse_scale = 1e4
    elif "lyapunov" in dataset.file_name:
        ellipse_scale = 5e5
    
    # Thrust
    thrust_color = "#e63946"
    thurst_alpha = 0.8
    show_thrust = True
    
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
    list_colors_reference = ["#7f7f7f", "#7f7f7f"]
    list_linestyles_references = ["dotted", "dashed"]
    
    # Trajectory
    color_trajectory = "black"
    linewidth_trajectory = 1.5

    # Sample
    color_sample = "#9a7bb5"
    max_sample_size = 200
    sample_alpha = 0.3
    sample_linewidth = 0.5
    
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

    # Get data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Nominal state"):
            data_state = dataset.list_datasets[i].copy()
            list_names_state = dataset.list_dataset_names[i].copy()

    for axies in list_axis:
        axis_0 = axies[0]
        axis_1 = axies[1]
        coord_0 = data_state[axis_0,:]
        coord_1 = data_state[axis_1,:]

        # Get sample
        sample_size = min(max_sample_size, int(len(dataset_sample.list_dataset_names)/4))
        coord_0_sample = []
        coord_1_sample = []
        if sample_size != 0:
            nb_datadets = len(dataset_sample.list_dataset_names)
            for i in range(0, sample_size):
                for j in range(nb_datadets):
                    if (dataset_sample.list_dataset_names[j][0] == "Sample state " + str(i)):
                        data_control_sample = dataset_sample.list_datasets[j].copy()
                coord_0_sample.append(data_control_sample[axis_0,:])
                coord_1_sample.append(data_control_sample[axis_1,:])

        # Normalisation
        if denormalise:
            LU = dataset.spacecraft_parameters.constants.LU
            coord_0 *= LU
            coord_1 *= LU

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
            
            # Sample + scaling
            if sample_size != 0:
                for i in range(sample_size):
                    coord_0_sample[i] = interpolate.interp1d(
                        t_old, coord_0 + ellipse_scale*(coord_0_sample[i] - coord_0), kind='cubic')(t_new)
                    coord_1_sample[i] = interpolate.interp1d(
                        t_old, coord_1 + ellipse_scale*(coord_1_sample[i] - coord_1), kind='cubic')(t_new)

            coord_0 = interpolate.interp1d(t_old, coord_0, kind='cubic')(t_new)
            coord_1 = interpolate.interp1d(t_old, coord_1, kind='cubic')(t_new)

        # Create plot
        fig = plt.figure(dpi=dpi)
        ax = fig.add_subplot()
        ax.set_aspect("equal")

        # Set labels
        ax.set_xlabel(list_names_state[axis_0 + 1])
        ax.set_ylabel(list_names_state[axis_1 + 1])

        # Sample
        if sample_size != 0:
            for i in range(sample_size):
                ax.plot(coord_0_sample[i], coord_1_sample[i],
                    alpha=sample_alpha,
                    linewidth=sample_linewidth,
                    color=color_sample,
                    zorder=1)
        
        # Plot state dispersion
        plot_state_distribution(dataset, axis_0, axis_1, ax,
            plot_CL, ellipse_scale, ellipse_alpha,
            ellipse_color, ellipse_linewidth, denormalise)
        
        # Plot reference orbits
        plot_reference_orbits(dataset, axis_0, axis_1, ax,
                              list_colors_reference,
                              list_linestyles_references,
                              denormalise)
        
        # Plot departure and arrival points
        plot_departure_arrival(dataset, axis_0, axis_1, ax,
                               list_colors_departure_arrival,
                               list_markers_departure_arrival,
                               denormalise)
        
        # Plot system points
        plot_system_points(dataset, axis_0, axis_1, ax,
                           list_colors_system_points,
                           list_markers_system_points,
                           list_plots_system_points,
                           denormalise)
        
        # Plot Thrust
        if show_thrust:
            plot_thrust_vector(dataset, axis_0, axis_1,
                               ax, thurst_alpha, thrust_color,
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
            signature = ("_robust_2d_" + str(axis_0) + "_" + str(axis_1)
                + "." + saving_format)
            file_name = file_name.replace(
                ".dat", signature)
            
            # Save
            plt.savefig(file_name, bbox_inches='tight')    
           
        if show_plot:   
            plt.show()
        else:
            plt.close(fig)
    