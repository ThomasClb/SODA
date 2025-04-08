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

    index_t =0
    
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

    LAM_p = 0
    mean = np.array(list_center[k][:6,0])
    Sigma = list_data_Sigma[0][:,0].reshape((d,d))[:6,:6]
    sigma_tilde = 0.6715664864669252
    Sigma[3,3] = Sigma[3,3]/(sigma_tilde*sigma_tilde)
    Sigma[4,4] = Sigma[4,4]/(sigma_tilde*sigma_tilde)
    p = multivariate_normal(mean, Sigma)
    for i in range(sample_size):
        sample_i = list_sample[i]
        LAM_p = LAM_p + p.pdf(sample_i[:6,index_t])
    LAM_p = LAM_p/N

    print(LAM/LAM_p)

    return 