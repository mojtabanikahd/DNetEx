print("Library imported")

import random
import numpy as np

print("Library imported")

def get_delta_hat(cov_mat_a, cov_mat_b):
    precision_mat_a = np.linalg.inv(cov_mat_a)
    precision_mat_b = np.linalg.inv(cov_mat_b)
    return precision_mat_a - precision_mat_b


def format_float(x, d):
    return float(f"{x:.{d}f}")

def round_matrix(mat, d):
    vectorized_format_float = np.vectorize(format_float)
    formatted_matrix = vectorized_format_float(mat, d)
    return formatted_matrix

def get_kronecker_value(mat1, mat2, i, j):
    d = mat1.shape[0]
    return mat1[i//d, j//d] * mat2[i%d, j%d]

def get_kronecker_col(mat1, mat2, col):
    d = mat1.shape[0]
    return [get_kronecker_value(mat1, mat2, i, col) for i in range(d**2)]

def get_kronecker_multi_col(mat1, mat2, indices):
    col_list = []
    for idx in indices:
        col_list.append(get_kronecker_col(mat1, mat2, idx))
    col_sub_kron = np.array(col_list).T
    return col_sub_kron

def get_kronecker_sub_mat(all_rows_kronecker, indices):
    return all_rows_kronecker[np.ix_(indices)]

def get_lambda_hit(indices, cols_kron, val1, val2, v, lambda_t):
    final_result = np.zeros((cols_kron.shape[0], 2))

    cols_kron_idx = cols_kron[indices]  # Select rows corresponding to indices
    numerator = np.dot(cols_kron_idx, val1) - v[indices]

    denom1 = (1 - np.dot(cols_kron_idx, val2)).reshape(-1, 1)
    denom2 = (-1 - np.dot(cols_kron_idx, val2)).reshape(-1, 1)
    result1 = numerator / denom1
    result2 = numerator / denom2
    cat_results = np.concatenate((result2, result1), axis=1)
    final_result[indices] = cat_results
    cat_results = round_matrix(final_result, 5)
    cat_results[cat_results >= lambda_t] = -1

    max_value = np.max(cat_results)
    max_indices = np.where(cat_results == max_value)
    indices_sign = np.sign(max_indices[1] - .5)
    max_indices = list(zip(max_indices[0], indices_sign))
    sorted_indices = sorted(max_indices, key=lambda x: x[0])
    return max_value, sorted_indices

def get_lambda_cross(val1, val2, lambda_t):
    cross_vals = - val1/val2
    cross_vals = round_matrix(cross_vals, 5)
    cross_vals[cross_vals >= lambda_t] = -1
    max_value = np.max(cross_vals)
    max_indices = np.where(cross_vals == max_value)[0]
    return max_value, list(max_indices)

def construct_delta_hat(val1, val2, lam_val):
    return -val1 - lam_val * val2

def get_delta_hat_l1(val1, val2, lam_val):
    delta_hat = construct_delta_hat(val1, val2, lam_val)
    return np.abs(delta_hat).sum()

def alg_1(sigma, sigma_hat, active_set_pack, lambda_t):
    active_set = [i[0] for i in active_set_pack]
    sa = [i[1] for i in active_set_pack]
    d = sigma.shape[0]
    active_complement = [i for i in range(d**2) if i not in active_set]
    v = (sigma - sigma_hat).flatten('F').reshape(-1, 1)
    zero_mat = np.zeros_like(v)
    if len(active_set) == 0:
        v_val_cat = np.concatenate((v, -v), axis=1)
        v_val_cat = round_matrix(v_val_cat, 5)
        next_lambda = np.max(v_val_cat)
        max_indices = np.where(v_val_cat == next_lambda)
        indices_sign = np.sign(max_indices[1] - .5)
        max_indices = list(zip(max_indices[0], indices_sign))
        sorted_indices = sorted(max_indices, key=lambda x: x[0])
        return next_lambda, sorted_indices, 0
    else:
        va = v[np.ix_(active_set)]
        cols_kron_1 = get_kronecker_multi_col(sigma, sigma_hat, active_set)
        cols_kron_2 = get_kronecker_multi_col(sigma_hat, sigma, active_set)

        sub_kron_1 = get_kronecker_sub_mat(cols_kron_1, active_set)
        sub_kron_2 = get_kronecker_sub_mat(cols_kron_2, active_set)

        cols_kron = (cols_kron_1+cols_kron_2)/2
        sub_kron = (sub_kron_1 + sub_kron_2)/2

        sub_kron_inverse = np.linalg.inv(sub_kron)
        val1 = sub_kron_inverse @ va
        val2 = (sub_kron_inverse @ sa).reshape(-1, 1)

        lambda_hit_val, lambda_hit_idx = get_lambda_hit(active_complement, cols_kron, val1, val2, v, lambda_t)
        lambda_cross_val, lambda_cross_idx = get_lambda_cross(val1, val2, lambda_t)
        lambda_cross_idx = [active_set[i] for i in lambda_cross_idx]

        next_lambda = max(lambda_hit_val, lambda_cross_val)
        delta_hat_l1 = get_delta_hat_l1(val1, val2, next_lambda)
        delta_hat = construct_delta_hat(val1, val2, next_lambda)
        zero_mat[np.ix_(active_set)] = delta_hat
        delta_hat_mat = zero_mat.reshape(sigma.shape, order='F')

    if lambda_hit_val > lambda_cross_val:
        new_active_pack = sorted(active_set_pack + lambda_hit_idx, key=lambda x: x[0])
    else:
        new_active_set = [i for i in active_set if i not in lambda_cross_idx]
        new_active_pack = [i for i in active_set_pack if i[0] in new_active_set]

    return next_lambda, new_active_pack, delta_hat_mat

def construct_graph(dim, indices):
    graph = np.zeros((dim, dim))
    for idx in indices:
        graph[idx[0]//dim, idx[0]%dim] = idx[1]
    return graph

def generate_mirror_statistics_sgn(delta_hat_ols, delta_hat_lasso):
    ols_diag = np.diag(delta_hat_ols)
    empty_mat = np.zeros(delta_hat_lasso.shape)
    for i in range(delta_hat_ols.shape[0]):
        for i1 in range(delta_hat_ols.shape[1]):
            empty_mat[i, i1] = 2*delta_hat_ols[i, i1] + ols_diag[i] + ols_diag[i1]
    sign_mat = np.sign(delta_hat_lasso * empty_mat)
    return sign_mat

def get_e_value(mirror_statistics, cutoff_value, dimension):
    sub_mat = list(mirror_statistics[np.triu_indices(dimension,k=0)])
    t1 = len(sub_mat)
    t2 = (mirror_statistics >= cutoff_value).astype(int)
    t3 = 1 + len([i for i in sub_mat if i <= -cutoff_value])
    return t1 * t2 / t3

print("Library imported")

def run_lasso(lasso_s_cov_1, lasso_s_cov_2, ols_delta_hat, c=30):
    dimension = lasso_s_cov_1.shape[0]
    a_set_sa, lam = [], 1000000
    evals = []
    e_values = []
    delta_hats = []
    while len(a_set_sa) < c:
        new_lam, new_a_set_sa, delta_hat = alg_1(lasso_s_cov_1, lasso_s_cov_2, a_set_sa, lam)

        estimated_graph = construct_graph(dimension, a_set_sa)
        mirror_stats = generate_mirror_statistics_sgn(ols_delta_hat, estimated_graph)
        eval = np.max(get_e_value(mirror_stats, 1, dimension))
        evals.append(eval)
        e_values.append(get_e_value(mirror_stats, 1, dimension))
        delta_hats.append(delta_hat)

        a_set_sa, lam = new_a_set_sa, new_lam

        if new_lam <= 0:
            break
        
    return delta_hats, e_values, evals