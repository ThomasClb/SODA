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

def compute_LAM(dataset, dataset_sample):

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_alpha = []
    list_data_Sigma = []
    list_center = []
    quad = np.zeros((2,2))
    x_min, x_max = 1e15, -1e15
    y_min, y_max = 1e15, -1e15
    for k in range(nb_GMM):
        for i in range(nb_datasets):
            if (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
                list_center.append(data_state)
            elif (dataset.list_dataset_names[i][0] == "Sigma GMM " + str(k)):
                data_Sigma = dataset.list_datasets[i].copy()
                d = int(np.sqrt(len(data_Sigma[:,0])))
                N = len(data_Sigma[0,:])
                list_data_Sigma.append(data_Sigma)
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

    index_t = N - 1
    
    list_pg = []
    for k in range(nb_GMM):
        centers = list_center[k]
        Sigmas = list_data_Sigma[k]
        pg_k = multivariate_normal(centers[:6,index_t], Sigmas[:,index_t].reshape((d,d))[:6,:6])
        list_pg.append(pg_k)

    # Get sample
    sample_size = int(len(dataset_sample.list_dataset_names)/4)
    list_sample = []
    if sample_size != 0:
        nb_datasets = len(dataset_sample.list_dataset_names)
        for i in range(0, sample_size):
            for j in range(nb_datasets):
                if (dataset_sample.list_dataset_names[j][0] == "Sample state " + str(i)):
                    data_control_sample = dataset_sample.list_datasets[j].copy()
            list_sample.append(data_control_sample)

   
    # Sample
    LAM = 0
    if sample_size != 0:
        for k in range(nb_GMM):
            alpha_k = list_alpha[k]
            pg_k = list_pg[k]
            LAM_k = 0.0
            for i in range(sample_size):
                sample_i = list_sample[i]
                LAM_k = LAM_k + pg_k.pdf(sample_i[:6,index_t])
            LAM = LAM + alpha_k*LAM_k
    LAM = LAM/N
    print(LAM)

    return 

"""
    Plots a 2D plot from a given dataset.

"""
def plot_2d_pdf(dataset, dataset_sample=Dataset()):
    
    # Settings
    plt.rcParams.update({'font.size': 10})
    dpi = 200
    
    # Axes
    list_axis = [[0, 1], [3, 4]]
    if "halo" in dataset.file_name:
        list_axis = [[0, 1], [0, 2]]
    sampling = 4
    
    # Ellipses
    cmap = "plasma" 
    ellipse_alpha = 1.0
    ellipse_nb_points = 604
    ellipse_levels = np.array([10**i for i in range(-5, 0 + 1, 1)])
    plot_CL = True
    if "mars" in dataset.file_name:
        ellipse_scale = 15
    elif "halo" in dataset.file_name:
        ellipse_scale = 3e5
    elif "nrho" in dataset.file_name:
        ellipse_scale = 2e5
    elif "dro_to_dro" in dataset.file_name:
        ellipse_scale = 1e4
    elif "lyapunov" in dataset.file_name:
        ellipse_scale = 5e5
    
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
    maker_sample = "."
    max_sample_size = 400
    sample_alpha = 0.8
    marker_size = 10
    
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
    show_plot = True

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
        fig, ax = plt.subplots(shape[0], shape[1], constrained_layout = True, dpi=dpi)
        for i, index in enumerate(list_index):
            ax_i = ax.flat[i]

            # Set labels
            ax_i.ticklabel_format(useMathText=True)
            
            x_min, x_max, y_min, y_max = plot_sample(
                dataset, dataset_sample, axis_0, axis_1, ax_i, index, 
                plot_CL, max_sample_size,
                maker_sample, sample_alpha, color_sample, marker_size,
                denormalise, interpolation, interpolation_rate)
            
            # Plot state dispersion TO DO
            CS = plot_state_distribution(dataset, axis_0, axis_1, ax_i, index,
                plot_CL, ellipse_nb_points,
                ellipse_alpha, ellipse_levels, [x_min, x_max, y_min, y_max],
                cmap, denormalise)

            # Plot departure and arrival points
            plot_departure_arrival(dataset, axis_0, axis_1, ax_i,
                                   index,
                                   list_colors_departure_arrival,
                                   list_markers_departure_arrival,
                                   denormalise)

            # lims
            x_m, y_m = 0.5*(x_min + x_max), 0.5*(y_min + y_max)
            dx, dy = (x_max - x_min)/2, (y_max - y_min)/2
            coef = 1.3
            ax_i.set_xlim(x_m - coef*dx, x_m + coef*dx)
            ax_i.set_ylim(y_m - coef*dy, y_m + coef*dy)

            step_str = r"$t/ToF=" + str(float("{:.2f}".format((1.0*index)/N))) + "$"
            ax_i.set_title(step_str)

            if show_grid:
                ax_i.grid(alpha=0.5, color="#d3d3d3")

        # Color bar
        cbar = fig.colorbar(CS, ax=ax, location="right")
        cbar.set_label("Normalized state PDF [-]")


        fig.text(0.5, 0.04, list_names_state[axis_0 + 1], ha='center', va='center')
        fig.text(0.06, 0.5, list_names_state[axis_1 + 1], ha='center', va='center', rotation='vertical')
        
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
    