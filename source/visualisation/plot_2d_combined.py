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
import os, sys
import subprocess


plt.rcParams.update({'font.size': 15})
ALPHA_0_GMM = 0.5495506294920584 # Central weight [-]
ALPHA_1_GMM = 0.225224685253970 # Lateral weight [-]



def plot_system_points_traj(dataset, axis_0, axis_1, ax,
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
def plot_departure_arrival_traj(dataset, axis_0, axis_1, ax,
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
        data_departure *= LU
        data_arrival *= LU
            
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
def plot_thrust_vector_traj(dataset, axis_0, axis_1, ax,
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
def plot_reference_orbits_traj(dataset, axis_0, axis_1, ax,
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
def plot_thrust_vector_traj(dataset, axis_0, axis_1, ax, index,
    alpha, color, denormalise):
    # Retreive data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Nominal state GMM 0"):
            data_state = dataset.list_datasets[i].copy()
        elif (dataset.list_dataset_names[i][0] == "Nominal control GMM 0"):
            data_control = dataset.list_datasets[i].copy()
    coord_0 = data_state[axis_0, index]
    coord_1 = data_state[axis_1, index]
    ucoord_0 = 10*data_control[axis_0,index]
    ucoord_1 = 10*data_control[axis_1,index]
    
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        coord_0 *= LU
        coord_1 *= LU
        ucoord_0 *= THRUSTU
        ucoord_1 *= THRUSTU
    
    # Plot arrows
    ax.quiver(coord_0, coord_1,
              ucoord_0, ucoord_1,
              color=color, label='Thrust',
              alpha=alpha, scale=5,
              zorder=10)
def plot_traj(dataset,
            ax, denormalise, show_grid, interpolation_rate,
            axis_0, axis_1, coord_0, coord_1, coord_0_inter, coord_1_inter,
            index, N):
    
    # Settings

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

    # Trajectory
    color_trajectory_future = "black"
    linestyle_trajectory_future = "dotted"
    color_trajectory = "black"
    linewidth_trajectory = 5
    linewidth_trajectory_future = 1.5
    marker_color = "#9A7BB5"
    marker = "d"
    marker_size = 75
    
    # Thrust
    thrust_color = "#e63946"
    thurst_alpha = 0.8
    show_thrust = True

    # Create plot
    ax.set_aspect("equal")

    # Plot reference orbits
    plot_reference_orbits_traj(dataset, axis_0, axis_1, ax,
                        alpha_reference,
                        list_colors_reference,
                        list_linestyles_references,
                        denormalise)
    
    # Plot departure and arrival points
    plot_departure_arrival_traj(dataset, axis_0, axis_1, ax,
                        list_colors_departure_arrival,
                        list_markers_departure_arrival,
                        denormalise)

    # Plot Thrust
    if show_thrust and index != N:
        plot_thrust_vector_traj(dataset, axis_0, axis_1,
                            ax, index, thurst_alpha, thrust_color,
                        denormalise)
    
    # Plot system points
    plot_system_points_traj(dataset, axis_0, axis_1, ax,
                    list_colors_system_points,
                    list_markers_system_points,
                    list_plots_system_points,
                    denormalise)

    # Plot trajectory
    scale = (N+1.0)/N
    ax.plot(coord_0_inter[:int(index*interpolation_rate*scale)], coord_1_inter[:int(index*interpolation_rate*scale)],
            color=color_trajectory,
            label='Trajectory',
            zorder=11,
            linewidth=linewidth_trajectory)
    ax.scatter(coord_0[index], coord_1[index],
            color=marker_color,
            marker=marker,
            s=marker_size,
            zorder=12)
    ax.plot(coord_0_inter[int(index*interpolation_rate*scale):], coord_1_inter[int(index*interpolation_rate*scale):],
            color=color_trajectory_future,
            zorder=11,
            linestyle=linestyle_trajectory_future,
            linewidth=linewidth_trajectory_future)

    if show_grid:
        ax.grid(alpha=0.5, color="#d3d3d3")



"""
	plot_2d.py

	Purpose: Implements the functions to produce 2D plots of tranfers.

	@author Thomas Caleb

	@version 2.0 12/12/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter
import scipy.interpolate as interpolate
from matplotlib.patches import Ellipse
from scipy.stats import chi2, multivariate_normal
plt.rcParams['hatch.linewidth'] = 0.5  # Default is ~1.0
import matplotlib.gridspec as gridspec

from misc import get_Lagrange_point
from classes import Dataset
from mpl_toolkits.axisartist.axislines import Subplot

ALPHA_0_GMM = 0.5495506294920584 # Central weight [-]
ALPHA_1_GMM = 0.225224685253970 # Lateral weight [-]

        
  
def plot_departure_arrival_pdf(dataset, axis_0, axis_1, ax, index,
                           list_colors, list_markers, lims,
                           transparancy,
                           denormalise):    
    # Retrieve data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Sigma GMM 0"):
                data_Sigma = dataset.list_datasets[i].copy()
                N = len(data_Sigma[0,:]) - 1
        if (dataset.list_dataset_names[i][0] == "Arrival Sigma"):
                data_arrival_Sigma = dataset.list_datasets[i].copy()
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i].copy()
        if (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i].copy()
    d = int(np.sqrt(len(data_Sigma[:,0])))
                
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        data_arrival *= LU
        data_departure *= LU
        data_arrival_Sigma *= LU*LU
            
    # Plot departure and arrival
    if index == 0:
        ax.scatter(data_departure[axis_0, 0], data_departure[axis_1, 0],
            s=10,
            color=list_colors[0], marker=list_markers[0],
                  zorder=100)
        ax.text(data_departure[axis_0, 0], data_departure[axis_1, 0],
                " $x_0$",
                  zorder=100)

    if index == N:
        if (data_arrival[axis_0, 0] > lims[0] and data_arrival[axis_0, 0] < lims[1] and
            data_arrival[axis_1, 0] > lims[2] and data_arrival[axis_1, 0] < lims[3]):
            ax.scatter(data_arrival[axis_0, 0], data_arrival[axis_1, 0],
                s=10,
                color=list_colors[1], marker=list_markers[1],
                      zorder=100)
            ax.text(data_arrival[axis_0, 0], data_arrival[axis_1, 0],
                    " $x_t$",
                      zorder=100)

        # Store covariance
        center = np.array([data_arrival[axis_0,0], data_arrival[axis_1,0]])
        Sigma = np.array(data_arrival_Sigma[:,0]).reshape((d-2,d-2))
        quad = np.zeros((2,2))
        quad[0,0] = Sigma[axis_0, axis_0]
        quad[1,1] = Sigma[axis_1, axis_1]
        quad[0,1] = Sigma[axis_0, axis_1]
        quad[1, 0] = quad[0,1]

        # Plot ellipses
        lambda_, v = np.linalg.eig(quad)
        lambda_ = 2*np.sqrt(chi2.ppf(0.95, d - 2)*lambda_)
        ell = Ellipse(
                 xy=center,
                  width=lambda_[0], height=lambda_[1],
                  angle=np.rad2deg(np.arccos(v[0, 0])),
                  alpha=transparancy,
                  hatch="///////////////",
                  fill=False,
                  linewidth=0.8,
                  color=list_colors[2],
                  zorder=-2)
        ax.add_artist(ell)
def plot_state_distribution_pdf(dataset, axis_0, axis_1, ax, index_, plot_CL, nb_points, coefficient_lims,
                            transparancy, levels, lims,
                            cmap, denormalise):

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_data_state = []
    list_data_Sigma = []
    list_data_der = []
    list_history = []
    for k in range(nb_GMM):
        for i in range(nb_datasets):
            if (dataset.list_dataset_names[i][0] == "Sigma GMM " + str(k)):
                data_Sigma = dataset.list_datasets[i].copy()
                list_data_Sigma.append(data_Sigma)
            elif (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
                list_data_state.append(data_state)
            elif (dataset.list_dataset_names[i][0] == "Der dynamics GMM " + str(k)):
                data_der = dataset.list_datasets[i].copy()
                list_data_der.append(data_der)
            elif (dataset.list_dataset_names[i][0] == "Splitting history GMM " + str(k)):
                data_history = dataset.list_datasets[i].copy()
                list_history.append(data_history)
    LU = dataset.spacecraft_parameters.constants.LU
        
    # Plot ellispses
    quad = np.zeros((2,2))
    d = int(np.sqrt(len(data_Sigma[:,0])))
    N = len(data_Sigma[0,:])

    i = index_
    Z_i = np.zeros((nb_points,nb_points))
    list_alpha = []
    list_pg = []
    # Get mean + cov
    for k in range(nb_GMM):            
        # Get alpha
        alpha = 1
        data_history = list_history[k]
        if data_history.shape[1] != 0:
            for j in range(len(data_history[0,:])):
                if data_history[1,j] == 0:
                    alpha = alpha*ALPHA_0_GMM
                elif abs(data_history[1,j]) == 1:
                    alpha = alpha*ALPHA_1_GMM
        list_alpha.append(alpha)
        center = np.array([list_data_state[k][:,i][axis_0], list_data_state[k][:,i][axis_1]])

        # Store covariance
        Sigma = list_data_Sigma[k][:,i]
        quad[0,0] = Sigma[axis_0 + d*axis_0]
        quad[1,1] = Sigma[axis_1 + d*axis_1]
        quad[0,1] = Sigma[axis_0 + d*axis_1]
        quad[1, 0] = quad[0,1]

        # Denormalise
        if denormalise:
            center *= LU
            quad *= LU*LU

        list_pg.append(multivariate_normal(mean=center, cov=quad))

    x_m, y_m = 0.5*(lims[0] + lims[1]), 0.5*(lims[2] + lims[3])
    dx, dy = (lims[1] - lims[0])/2, (lims[3] - lims[2])/2
    coef = coefficient_lims
    X_norm, Y_norm = np.meshgrid(np.linspace(x_m - coef*dx, x_m + coef*dx, nb_points), np.linspace(y_m - coef*dy, y_m + coef*dy, nb_points))
    pos_norm = np.dstack((X_norm, Y_norm))
    X = pos_norm[0,:,0]
    Y = pos_norm[:,0,1]

    # Denormalise
    if denormalise:
        X *= LU
        Y *= LU

    for k in range(nb_GMM):                    
        # Add pdfs
        Z_i = Z_i + list_alpha[k]*list_pg[k].pdf(pos_norm)
    Z_i /= np.max(Z_i)

    CS = ax.contourf(X, Y, Z_i, levels,
        alpha=transparancy,
        cmap=cmap, zorder=-1)

    return CS
def plot_sample_pdf(dataset, dataset_sample, axis_0, axis_1, ax, index,
                plot_CL, max_sample_size,
                maker, transparancy, color, marker_size, linewidth,
                denormalise, interpolation, interpolation_rate):

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_alpha = []
    list_data_Sigma = []
    list_inv_cov = []
    quad = np.zeros((2,2))
    x_min, x_max = 1e15, -1e15
    y_min, y_max = 1e15, -1e15
    for k in range(nb_GMM):
        for i in range(nb_datasets):
            if (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
                list_names_state = dataset.list_dataset_names[i].copy()
            elif (dataset.list_dataset_names[i][0] == "Sigma GMM " + str(k)):
                data_Sigma = dataset.list_datasets[i].copy()
                d = int(np.sqrt(len(data_Sigma[:,0])))
                N = len(data_Sigma[0,:])

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
        
    # Sample
    if sample_size != 0:
        for i in range(sample_size):
            x = coord_0_sample[i][index]
            y = coord_1_sample[i][index]
            if x > x_max:
                x_max = x
            if x < x_min:
                x_min = x
            if y > y_max:
                y_max = y
            if y < y_min:
                y_min = y
            ax.scatter(x, y,
                alpha=transparancy,
                s=marker_size,
                marker=maker,
                color=color,
                linewidth=linewidth,
                zorder=10)
    return x_min, x_max, y_min, y_max
def plot_pdf(dataset, dataset_sample, fig, ax_i, cax, 
             axis_0, axis_1, denormalise, interpolation, interpolation_rate, show_grid,
             index, N):
    
    # Settings
    coefficient_lims = 1.3
    
    # Ellipses
    cmap = "plasma" 
    ellipse_alpha = 0.8
    ellipse_nb_points = 404
    ellipse_levels =  [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0] # np.linspace(0.01, 1, nb_levels)
    plot_CL = True
    
    # Departure arrival
    alpha_departure_arrival = 0.45
    list_colors_departure_arrival = ["black", "black","#4caf50"]
    list_markers_departure_arrival = ["^", "v"]
    

    # Sample
    color_sample = "#e57373"
    color_sample = "white"
    maker_sample = "+"
    max_sample_size = 400
    sample_alpha = 0.7
    marker_size = 30
    marker_linewidth = 1

    # Locators
    ax_i.set_aspect("equal")
    ax_i.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax_i.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax_i.tick_params(axis='both', which='major') #, labelsize=10)

    # Formatter
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0,0))
    formatter.set_useOffset(True)
    ax_i.xaxis.set_major_formatter(formatter)
    ax_i.yaxis.set_major_formatter(formatter)

    # 🔑 Fixer une largeur minimale pour les labels d’offset
    ax_i.xaxis.get_offset_text().set_size(10)
    ax_i.yaxis.get_offset_text().set_size(10)
    ax_i.get_xaxis().get_offset_text().set_x(0.95)  # place fixe
    ax_i.get_yaxis().get_offset_text().set_y(1.02)  # place fixe

    # Réserver de l’espace pour ne pas que ça bouge
    ax_i.ticklabel_format(style='sci', axis='both', scilimits=(0,0), useMathText=True)

    
    lims = plot_sample_pdf(
        dataset, dataset_sample, axis_0, axis_1, ax_i, index, 
        plot_CL, max_sample_size,
        maker_sample, sample_alpha, color_sample, marker_size,
        marker_linewidth,
        denormalise, interpolation, interpolation_rate)
    
    # Plot state dispersion
    coef = coefficient_lims
    CS = plot_state_distribution_pdf(dataset, axis_0, axis_1, ax_i, index,
        plot_CL, ellipse_nb_points, coef,
        ellipse_alpha, ellipse_levels, lims,
        cmap, denormalise)

    # Plot departure and arrival points
    plot_departure_arrival_pdf(dataset, axis_0, axis_1, ax_i,
                            index, 
                            list_colors_departure_arrival,
                            list_markers_departure_arrival, lims,
                            alpha_departure_arrival,
                            denormalise)

    # lims
    x_min, x_max, y_min, y_max = lims
    x_m, y_m = 0.5*(x_min + x_max), 0.5*(y_min + y_max)
    dx, dy = (x_max - x_min)/2, (y_max - y_min)/2
    d_max = max(dx, dy)
    ax_i.set_xlim(x_m - coef*d_max, x_m + coef*d_max)
    ax_i.set_ylim(y_m - coef*d_max, y_m + coef*d_max)
    if (index == 0 or index == N):
        step_str = r"$t/ToF=" + str(int((1.0*index)/N)) + "$"
    else:
        step_str = r"$t/ToF=" + str(float("{:.2f}".format((1.0*index)/N))) + "$"
    ToF_text = ax_i.text(x_m - 0.8*coef*d_max, y_m + 0.8*coef*d_max,
            step_str,
            fontsize="xx-large",
            zorder=100)
    ToF_text.set_bbox(dict(
        facecolor=(1, 1, 1, 0.3),   # white with X% opacity
        edgecolor=(0, 0, 0, 1.0),   # black with full opacity
        linewidth=1.0,
        boxstyle='round,pad=0.3'  # Optional: rounded corners and padding
    ))

    if show_grid:
        ax_i.grid(alpha=0.5, color="#d3d3d3")

    # Color bar
    cbar = fig.colorbar(CS, cax=cax)    
    cbar.set_label("PDF [-]")
    

def plot_gif(dataset, dataset_sample=Dataset()):
    
    # Settings
    
    dpi = 400
    
    # Axes
    list_axis = [[0, 1]]
    if "halo" in dataset.file_name:
        list_axis = [[0, 2]]

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
    saving_format = "png"
    show_plot = False

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_data_state = []
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
        N = len(coord_0) - 1

        # Normalisation
        if denormalise:        
            # Labels
            list_names_state[axis_0 + 1] = list_names_state[axis_0 + 1].replace(
                "LU", "km")
            list_names_state[axis_1 + 1] = list_names_state[axis_1 + 1].replace(
                "LU", "km")

        if interpolation:
            t_old = np.linspace(0, 1, N + 1)  
            t_new = np.linspace(0, 1, interpolation_rate*(N + 1))  
            coord_0_inter = interpolate.interp1d(t_old, coord_0, kind='cubic')(t_new)
            coord_1_inter = interpolate.interp1d(t_old, coord_1, kind='cubic')(t_new)


        np_repeat_final = 3
        np_repeat_first = 3
        for index in range(0, np_repeat_first + N + 1 + np_repeat_final):

            # Create plot
            fig = plt.figure(figsize=(20, 10), constrained_layout=True, dpi=dpi)
            gs = gridspec.GridSpec(1, 3, width_ratios=[20, 20, 1], wspace=0.3, hspace=0.3)  # 20 parts plot, 1 part cbar
            ax_traj = fig.add_subplot(gs[0])
            ax_pdf = fig.add_subplot(gs[1])
            cax = fig.add_subplot(gs[2]) 

            if index >= N + np_repeat_first:
                index_ = N
            elif index < np_repeat_first:
                index_ = 0
            else:
                index_ = index - np_repeat_first

            plot_traj(dataset,
                ax_traj, denormalise, show_grid, interpolation_rate,
                axis_0, axis_1, coord_0, coord_1, coord_0_inter, coord_1_inter,
                index_, N)
            
            plot_pdf(dataset, dataset_sample, fig, ax_pdf, cax, 
                axis_0, axis_1, denormalise, interpolation, interpolation_rate, show_grid,
                index_, N)
            
            fig.text(0.45, 0.12, list_names_state[axis_0 + 1], ha='center', va='center', fontsize=15)
            fig.text(0.09, 0.5, list_names_state[axis_1 + 1], ha='center', va='center', rotation='vertical', fontsize=15)

            # Layout
            # fig.tight_layout(pad=0.2)
            
            if show_legend: 
                plt.legend(loc=legend_loc)
            
            if save_figure:
                # Get adress
                file_name = dataset.file_name.replace(
                    "robust_trajectory", "plots/gif")
                        
                # Add signature
                signature = "_" + (str(axis_0) + "_" + str(axis_1) + "/" + str(index)
                    + "." + saving_format)
                file_name = file_name.replace(
                    ".dat", signature)
                
                # Save
                plt.savefig(file_name, bbox_inches='tight')    
            
            if show_plot:   
                plt.show()
            else:
                plt.close(fig)
    # os.system("ffmpeg -start_number 0 -framerate 6 -i "cr3bp_EARTH_MOON_lt_haloL2_to_haloL1_20_20_5e-4_20_0_2_%d.png" -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -pix_fmt yuv420p output.mp4")  
