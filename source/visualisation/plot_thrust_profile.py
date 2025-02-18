"""
    plot_thrust_profile.py

    Purpose: Implements the functions to produce thrust profiles of tranfers.

    @author Thomas Caleb

    @version 2.0 12/12/2024
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from scipy.stats import chi2, multivariate_normal

from classes import Dataset

plt.rcParams.update({'font.size': 15})
ALPHA_0_GMM = 0.5495506294920584 # Central weight [-]
ALPHA_1_GMM = 0.225224685253970 # Lateral weight [-]

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


def plot_sample(dataset, dataset_sample, ax, max_sample_size,
                scale, transparancy, color, linewidth,
                denormalise):

    # Retreive data
    nb_datasets = len(dataset.list_dataset_names)

    # Normalisation
    if denormalise:
        THRUSTU = dataset.spacecraft_parameters.constants.THRUSTU
        TU = dataset.spacecraft_parameters.constants.TU

    # Retrieve GMM size
    nb_GMM = int((nb_datasets - 2)/6)
    list_control = []
    list_gains = []
    list_alpha = []
    list_data_Sigma = []
    list_center = []
    list_inv_cov = []
    quad = np.zeros((6,6))
    for k in range(nb_GMM):
        for i in range(nb_datasets):
            if (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
                list_names_state = dataset.list_dataset_names[i].copy()
                list_center.append(data_state[:6,0])
            elif (dataset.list_dataset_names[i][0] == "Nominal control GMM " + str(k)):
                data_control = dataset.list_datasets[i].copy()
                ux = data_control[0,:]
                uy = data_control[1,:]
                uz = data_control[2,:]
                dt = data_state[7,:-1]
                N = len(data_state[7,:])
                u = np.sqrt(ux*ux + uy*uy + uz*uz)
                list_control.append(u)
            elif (dataset.list_dataset_names[i][0] == "Feedback gain GMM " + str(k)):
                data_gains = dataset.list_datasets[i].copy()
                list_gains.append(data_gains)
            elif (dataset.list_dataset_names[i][0] == "Sigma GMM " + str(k)):
                data_Sigma = dataset.list_datasets[i].copy()
                d = int(np.sqrt(len(data_Sigma[:,0])))
                N = len(data_Sigma[0,:])
                list_data_Sigma.append(data_Sigma)
                Sigma = data_Sigma[:,0]
                for m in range(6):
                    quad[m,m] = Sigma[m + d*m]
                    for n in range(m):
                        quad[m,n] = Sigma[m + d*n]
                        quad[m,n] = quad[n,m]
                list_inv_cov.append(np.linalg.inv(quad))

    # Get sample
    sample_size = min(max_sample_size, int(len(dataset_sample.list_dataset_names)/4))
    control_sample = []
    state_0_sample = []
    if sample_size != 0:
        nb_datasets = len(dataset_sample.list_dataset_names)
        for i in range(0, sample_size):
            for j in range(nb_datasets):
                if (dataset_sample.list_dataset_names[j][0] == "Sample control " + str(i)):
                    data_control_sample = dataset_sample.list_datasets[j].copy()
                    ux = data_control_sample[0,:]
                    uy = data_control_sample[1,:]
                    uz = data_control_sample[2,:]
                    u = np.sqrt(ux*ux + uy*uy + uz*uz)
                    control_sample.append(u)
                if (dataset_sample.list_dataset_names[j][0] == "Sample state " + str(i)):
                    data_state_sample = dataset_sample.list_datasets[j].copy()
                    state_0_sample.append(data_state_sample[:6,0])

    # Normalisation
    if denormalise:
        LU = dataset.spacecraft_parameters.constants.LU
        # Normalisation
        if denormalise:
            for k in range(nb_GMM):
                list_control[k] = list_control[k]*THRUSTU
            dt *= (TU/86400)

        # Sample
        if sample_size != 0:
            for i in range(sample_size):
                control_sample[i] = THRUSTU*control_sample[i]
        
    # Offset
    if sample_size != 0:
        for i in range(sample_size):
            # Unpack
            u = control_sample[i]
            state_0 = state_0_sample[i]

            # Find split
            list_maha = []
            for k in range(nb_GMM):
                list_maha.append((state_0 - list_center[k]).T@list_inv_cov[k]@(state_0 - list_center[k]))
            index = np.argmin(list_maha)
            u_mean = list_control[index]

            # Scale
            control_sample[i] = u_mean + (u - u_mean)*scale**2

    # Sample
    if sample_size != 0:
        for i in range(sample_size):
            t, u_stairs = make_stairs(dt, control_sample[i])
            ax.plot(t, u_stairs,
                    color=color, linewidth=linewidth,
                    alpha=transparancy, zorder=1)


"""
    Plots a thrust pofile for a given transfer dataset.

"""
def plot_thrust_profile(dataset, dataset_sample=Dataset()):
    # Settings
    dpi = 200
    
    # Thrust norm
    color_thrust_norm = "#9a7bb5"
    label_thrust_norm  = "Thrust"
    alpha_thrust_norm = 0.5
    linewidth_thrust_norm = 1
    sampling_thrust = 10000
     
    # Max thrust norm
    color_thrust_max = "#e63946"
    linestyle_thrust_max = "dotted"
    label_thrust_max  = "Max thrust"
    alpha_thrust_max = 0.5

    # Robust thrust norm
    color_thrust_robust = "#4a90e2"
    linestyle_thrust_robust = "dashed"
    alpha_thrust_robust = 0.5
    linewidth_thrust_robust = 0.8
    label_thrust_robust  = "Margins"
    if "mars" in dataset.file_name:
        scale_thrust_robust = 10
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

    # Define a colormap
    cmap = "plasma_r" 
    
    # Normalisation
    denormalise = True
    
    # Legend
    show_legend = False
    legend_loc = "lower right"
    x_label = "Time [TU]"
    y_label = "Thrust norm [THRUSTU]"

    # Output
    show_grid = True
    save_figure = True
    saving_format = "pdf"
    show_plot = False
    
    # Retreive data
    nb_datadets = len(dataset.list_dataset_names)
    max_trust = dataset.spacecraft_parameters.thrust

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
        max_trust *= THRUSTU

        # Labels
        x_label = x_label.replace("TU", "days")
        y_label = y_label.replace("THRUSTU", "N")

    # Retrieve GMM size
    nb_GMM = int((nb_datadets - 2)/6)
    list_u = []
    x = np.linspace(0, 1.05*max_trust, sampling_thrust)
    list_pdf = []
    data_state_0 = []
    for k in range(nb_GMM):
        for i in range(nb_datadets):
            if (dataset.list_dataset_names[i][0] == "Sigma GMM " + str(k)):
                data_Sigma = dataset.list_datasets[i].copy()
            elif (dataset.list_dataset_names[i][0] == "Nominal control GMM " + str(k)):
                data_control = dataset.list_datasets[i].copy()
            elif (dataset.list_dataset_names[i][0] == "Nominal state GMM " + str(k)):
                data_state = dataset.list_datasets[i].copy()
            elif (dataset.list_dataset_names[i][0] == "Feedback gain GMM " + str(k)):
                data_gains = dataset.list_datasets[i].copy()
            elif (dataset.list_dataset_names[i][0] == "Splitting history GMM " + str(k)):
                data_history = dataset.list_datasets[i].copy()

        # Get mean
        ux = data_control[0,:]
        uy = data_control[1,:]
        uz = data_control[2,:]
        dt = data_state[7,:-1]
        N = len(data_state[7,:])
        u = np.sqrt(ux*ux + uy*uy + uz*uz)

        # Normalisation
        if denormalise:
            u *= THRUSTU
            ux *= THRUSTU
            uy *= THRUSTU
            uz *= THRUSTU
            dt *= (TU/86400)


        if k == 0:
            data_state_0 = data_state
            ux_0 = data_control[0,:]
            uy_0 = data_control[1,:]
            uz_0 = data_control[2,:]
            u_0 = np.sqrt(ux_0*ux_0 + uy_0*uy_0 + uz_0*uz_0)

        # Get Sigma
        list_sigma = []
        for j in range(N-1):
            K = data_gains[:,j].reshape((3,8))
            Sigma_x = data_Sigma[:,j].reshape((8,8))
            A = np.array([ux[j], uy[j], uz[j]]).T@K
            if u[j] != 0:
                A = A/u[j]
            Sigma = A@Sigma_x@A.T
            
            # Normalisation
            if denormalise:
                Sigma *= THRUSTU*THRUSTU
            
            # Scale
            Sigma = scale_thrust_robust*scale_thrust_robust*Sigma

            list_sigma.append(Sigma)

        # Get alpha
        alpha = 1
        if data_history.shape[1] != 0:
            for j in range(len(data_history[0,:])):
                if data_history[1,j] == 0:
                    alpha = alpha*ALPHA_0_GMM
                elif abs(data_history[1,j]) == 1:
                    alpha = alpha*ALPHA_1_GMM

        # Make pdfs
        l_pdf = []
        for j in range(N-1):
            l_pdf.append(alpha*multivariate_normal.pdf(x, mean=u[j], cov=list_sigma[j], allow_singular=False))
        list_pdf.append(l_pdf)

    thrust_pdf = np.zeros((sampling_thrust,N-1))
    for k in range(nb_GMM):
        for j in range(N-1):
            thrust_pdf[:,j] = thrust_pdf[:,j] + list_pdf[k][j]
    for j in range(N-1):
        thrust_pdf[:,j] /= max(thrust_pdf[:,j])
    thrust_pdf[thrust_pdf<1e-290] = 0

    # Create plot
    fig = plt.figure(dpi=dpi)
    ax = fig.add_subplot()

    # Set labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    # Axies limites
    plt.ylim((0 - 0.05*max_trust, 1.05*max_trust))
    plt.xlim((-dt[0], N*dt[0]))

    plot_sample(dataset, dataset_sample, ax, max_sample_size,
                scale_thrust_robust, sample_alpha, sample_color, sample_linewidth,
                denormalise)

    # Colormap
    im = ax.imshow(thrust_pdf,
        norm=clr.LogNorm(vmax=1,vmin=1e-250),
        aspect="auto",
        interpolation="None",
        origin="lower",
        alpha=alpha_thrust_robust,
        extent=[0, (N-1)*dt[0], x[0], x[-1]],
        cmap=cmap, zorder=10)
       
    # Plot max thrust
    t, u_stairs = make_stairs(dt, u*0 + max_trust)
    ax.plot(t, u_stairs, label=label_thrust_max,
            color=color_thrust_max, linestyle=linestyle_thrust_max,
            alpha=alpha_thrust_max)
    
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Normalized thrust norm PDF [-]")

    # Layout
    fig.tight_layout(pad=0.2)
    
    if show_grid:
        plt.grid(alpha=0.5, color="#d3d3d3",zorder=-1000)
    
    if show_legend: 
        plt.legend(loc=legend_loc)
    
    if save_figure:
        # Get adress
        file_name = dataset.file_name.replace(
            "robust_trajectory", "plots")
                
        # Add signature
        signature = ("_thrust" + "." + saving_format)
        file_name = file_name.replace(
            ".dat", signature)
        
        # Save
        plt.savefig(file_name, bbox_inches='tight')    
       
    if show_plot:   
        plt.show()
    else:
        plt.close(fig)
    