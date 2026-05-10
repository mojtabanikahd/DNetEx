import numpy as np
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri
# from rpy2.robjects.vectors import IntVector
from rpy2.robjects.vectors import FloatVector


# Output sequence is samples_A, precision_matrix_A, covariance_matrix_A,
# samples_B,  precision_matrix_B, covariance_matrix_B
def generate_reference_models(number_of_nodes, number_of_samples, number_of_changes, type="ScaleFree",
                              density_of_graph=0.2, power=1, mult=1):
    pandas2ri.activate()
    robjects.r.source('R_codes/Libraray.R')
    r_func = robjects.globalenv['generate_reference_models']
    output = r_func(number_of_nodes, number_of_samples, number_of_changes, type, density_of_graph, power, mult)
    samples_A, precision_mat_A, cov_mat_A = np.array(output[0]), np.array(output[1]), np.array(output[2])
    samples_B, precision_mat_B, cov_mat_B = np.array(output[3]), np.array(output[4]), np.array(output[5])
    return samples_A, precision_mat_A, cov_mat_A, samples_B, precision_mat_B, cov_mat_B

def DNetFinder_Liu2017(XA, XB, alphas, delta_star):
    pandas2ri.activate()
    numpy2ri.activate()
    # Convert Python lists or arrays to numpy arrays and ensure they are C-contiguous
    XA = np.ascontiguousarray(XA)
    XB = np.ascontiguousarray(XB)
    alphas = np.ascontiguousarray(alphas)
    delta_star = np.ascontiguousarray(delta_star)
  
    # Source your R code
    robjects.r.source('R_codes/Libraray.R')
    r_func = robjects.globalenv['DNetFinder_Liu2017']

    # Convert Python arrays to R matrices
    r_XA = robjects.r.matrix(XA, nrow=XA.shape[0], ncol=XA.shape[1])
    r_XB = robjects.r.matrix(XB, nrow=XB.shape[0], ncol=XB.shape[1])
    r_alphas = FloatVector(alphas)
    r_delta_star = robjects.r.matrix(delta_star, nrow=delta_star.shape[0], ncol=delta_star.shape[1])
  
    # Call the R function
    output = r_func(r_XA, r_XB, r_alphas, r_delta_star)
    output = pandas2ri.rpy2py(output)

    return output
