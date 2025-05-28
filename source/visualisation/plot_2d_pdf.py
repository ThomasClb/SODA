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
import scipy.interpolate as interpolate
from scipy.stats import chi2, multivariate_normal

from misc import get_Lagrange_point
from classes import Dataset


ALPHA_0_GMM = 0.5495506294920584 # Central weight [-]
ALPHA_1_GMM = 0.225224685253970 # Lateral weight [-]

        
"""
    Plots the departure and arrival points of a transfer.

"""    
def plot_departure_arrival(dataset, axis_0, axis_1, ax, index,
                           list_colors, list_markers, lims,
                           denormalise):    
    # Retrieve data
    nb_dataets = len(dataset.list_dataset_names)
    for i in range(nb_dataets):
        if (dataset.list_dataset_names[i][0] == "Sigma GMM 0"):
                data_Sigma = dataset.list_datasets[i].copy()
                N = len(data_Sigma[0,:]) - 1
        if (dataset.list_dataset_names[i][0] == "Departure orbit"):
            data_departure = dataset.list_datasets[i].copy()
        if (dataset.list_dataset_names[i][0] == "Arrival orbit"):
            data_arrival = dataset.list_datasets[i].copy()
                
    # Denormalise
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        data_state *= LU
            
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

"""
    Plots the projections of the uncertainty ellipsoïds of a transfer.

"""
def plot_state_distribution(dataset, axis_0, axis_1, ax, index_,
                            plot_CL, nb_points,
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
    coef = 1.3
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


def plot_sample(dataset, dataset_sample, axis_0, axis_1, ax, index,
                plot_CL, max_sample_size,
                maker, transparancy, color, marker_size,
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
                zorder=10)
    return x_min, x_max, y_min, y_max

"""
    Plots a 2D plot from a given dataset.

"""
def plot_2d_pdf(dataset, dataset_sample=Dataset()):
    
    # Settings
    plt.rcParams.update({'font.size': 8})
    dpi = 200
    
    # Axes
    list_axis = [[0, 1],[3, 4]]
    if "halo" in dataset.file_name:
        list_axis = [[0, 1], [0, 2]]
    sampling = 4
    
    # Ellipses
    cmap = "plasma" 
    ellipse_alpha = 0.8
    ellipse_nb_points = 404
    levels_min, levels_max = -5, 0
    nb_levels = 10
    ellipse_levels = np.linspace(0.1, 1, nb_levels)
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
    maker_sample = "+"
    max_sample_size = 400
    sample_alpha = 0.2
    marker_size = 30
    
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
                N = len(data_state[0,:]) - 1

    for axies in list_axis:
        axis_0 = axies[0]
        axis_1 = axies[1]

        # Normalisation
        if denormalise:        
            # Labels
            list_names_state[axis_0 + 1] = list_names_state[axis_0 + 1].replace(
                "LU", "km")
            list_names_state[axis_1 + 1] = list_names_state[axis_1 + 1].replace(
                "LU", "km")

        # Create plot
        list_index_plot = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
        list_index = [int(i/(sampling - 1)*N) for i in range(sampling)]
        shape = (2, 2)
        fig, ax = plt.subplots(shape[0], shape[1], dpi=dpi)
        fig.subplots_adjust(wspace=0.25, hspace=0.25)
        for i, index in enumerate(list_index):
            ax_i = ax.flat[i]

            # Set labels
            ax_i.ticklabel_format(useMathText=True)
            
            lims = plot_sample(
                dataset, dataset_sample, axis_0, axis_1, ax_i, index, 
                plot_CL, max_sample_size,
                maker_sample, sample_alpha, color_sample, marker_size,
                denormalise, interpolation, interpolation_rate)
            
            # Plot state dispersion TO DO
            CS = plot_state_distribution(dataset, axis_0, axis_1, ax_i, index,
                plot_CL, ellipse_nb_points,
                ellipse_alpha, ellipse_levels, lims,
                cmap, denormalise)

            # Plot departure and arrival points
            plot_departure_arrival(dataset, axis_0, axis_1, ax_i,
                                   index,
                                   list_colors_departure_arrival,
                                   list_markers_departure_arrival, lims,
                                   denormalise)

            # lims
            x_min, x_max, y_min, y_max = lims
            x_m, y_m = 0.5*(x_min + x_max), 0.5*(y_min + y_max)
            dx, dy = (x_max - x_min)/2, (y_max - y_min)/2
            coef = 1.3
            ax_i.set_xlim(x_m - coef*dx, x_m + coef*dx)
            ax_i.set_ylim(y_m - coef*dy, y_m + coef*dy)
            step_str = r"$t/ToF=" + str(float("{:.2f}".format((1.0*index)/N))) + "$"
            ax_i.text(x_m - dx, y_m + dy,
                    step_str,
                    zorder=100)

            if show_grid:
                ax_i.grid(alpha=0.5, color="#d3d3d3")

        # Color bar
        
        cbar = fig.colorbar(CS, ax=ax, location="right")
        cbar.set_label("PDF [-]")


        fig.text(0.5, 0.04, list_names_state[axis_0 + 1], ha='center', va='center')
        fig.text(0.06, 0.5, list_names_state[axis_1 + 1], ha='center', va='center', rotation='vertical')
        

        if show_legend:   
            plt.legend(loc=legend_loc)
        
        if save_figure:
            # Get adress
            file_name = dataset.file_name.replace(
                "robust_trajectory", "plots")
                    
            # Add signature
            signature = ("_2d_pdf_" + str(axis_0) + "_" + str(axis_1)
                + "." + saving_format)
            file_name = file_name.replace(
                ".dat", signature)
            
            # Save
            plt.savefig(file_name, bbox_inches='tight')    
           
        if show_plot:   
            plt.show()
        else:
            plt.close(fig)
    