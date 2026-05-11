from r_wrappers import SPDtrace, generate_reference_models, DNetFinder_Liu2017, DiffNetFDR_Xia2015, DiffNetFDR_Liu2017
import numpy as np
import random
import pickle
import time
import pandas as pd
import os
os.chdir(os.path.dirname(__file__))


np.random.seed(17)
random.seed(17)


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


def format_float(x, d):
    return float(f"{x:.{d}f}")


def round_matrix(mat, d):
    vectorized_format_float = np.vectorize(format_float)
    formatted_matrix = vectorized_format_float(mat, d)
    return formatted_matrix


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


def construct_graph(dim, indices):
    graph = np.zeros((dim, dim))
    for idx in indices:
        graph[idx[0]//dim, idx[0]%dim] = idx[1]
    return graph


def filter_unique_pairs(pairs):
    seen = set()
    unique_pairs = []

    for a, b in pairs:
        # Ensure the pair is always stored in a consistent order
        if (b, a) not in seen:
            seen.add((a, b))
            unique_pairs.append((a, b))

    return unique_pairs


def get_diff_nodes(mat):
    diff = np.where(mat != 0)
    diff_nodes = [(i, j) for i, j in zip(diff[0], diff[1])]
    diff_nodes = filter_unique_pairs(diff_nodes)
    flattened_list = [element for pair in diff_nodes for element in pair]
    return list(np.unique(flattened_list))


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


def split_dataset_samples(data_1, data_2, split_ratio=0.5):
    ols_samples_1 = random.sample(range(data_1.shape[0]), round(split_ratio*data_1.shape[0]))
    init_lasso_samples_1 = np.array(range(data_1.shape[0]))
    lasso_samples_1 = np.delete(init_lasso_samples_1, ols_samples_1, axis=0)
    ols_data_1 = data_1[np.ix_(ols_samples_1, )]
    lasso_data_1 = data_1[np.ix_(lasso_samples_1, )]

    ols_samples_2 = random.sample(range(data_2.shape[0]), round(split_ratio*data_2.shape[0]))
    init_lasso_samples_2 = np.array(range(data_2.shape[0]))
    lasso_samples_2 = np.delete(init_lasso_samples_2, ols_samples_2, axis=0)
    ols_data_2 = data_2[np.ix_(ols_samples_2, )]
    lasso_data_2 = data_2[np.ix_(lasso_samples_2, )]
    return lasso_data_1, lasso_data_2, ols_data_1, ols_data_2


def get_delta_hat(dataset_a, dataset_b):
    cov_mat_a = np.cov(dataset_a.T)
    precision_mat_a = np.linalg.inv(cov_mat_a)
    cov_mat_b = np.cov(dataset_b.T)
    precision_mat_b = np.linalg.inv(cov_mat_b)
    return precision_mat_a - precision_mat_b


def generate_mirror_statistics_sgn(delta_hat_ols, delta_hat_lasso):
    ols_diag = np.diag(delta_hat_ols)
    empty_mat = np.zeros(delta_hat_lasso.shape)
    for i in range(delta_hat_ols.shape[0]):
        for i1 in range(delta_hat_ols.shape[1]):
            empty_mat[i, i1] = 2*delta_hat_ols[i, i1] + ols_diag[i] + ols_diag[i1]
    sign_mat = np.sign(delta_hat_lasso * empty_mat)
    return sign_mat


def get_H0_H1(prc1, prc2):
    d_star = prc1 - prc2
    s_set_edges = np.argwhere(d_star != 0)
    H1_edges = [[i0, i1] for i0, i1 in s_set_edges.tolist() if i0 != i1]
    s_comp_edges = np.argwhere(d_star == 0)
    H0_edges = [[i0, i1] for i0, i1 in s_comp_edges.tolist() if i0 != i1]
    return H0_edges, H1_edges


def intersect_lists_of_lists(list1, list2):
    set1 = set(map(tuple, list1))
    set2 = set(map(tuple, list2))
    intersected_set = set1.intersection(set2)
    return [list(item) for item in intersected_set]


def get_fpr_tpr(H0_edges, H1_edges, mirror_statistics, cutoff_value):
    mirror_mat = mirror_statistics.copy()
    np.fill_diagonal(mirror_mat, 0)
    differential_elements = np.argwhere(mirror_mat >= cutoff_value)
    s_hat_cutoff = differential_elements.tolist()

    fp = intersect_lists_of_lists(s_hat_cutoff, H0_edges)
    tp = intersect_lists_of_lists(s_hat_cutoff, H1_edges)
    fn = [i for i in H1_edges if i not in s_hat_cutoff]

    if (len(tp) + len(fp)) == 0:
        precision = None
    else:
        precision = len(tp)/(len(tp) + len(fp))

    if (len(tp) + len(fn)) == 0:
        recall = None
    else:
        recall = len(tp)/(len(tp) + len(fn))

    if (len(fp) + len(tp)) == 0:
        fdr = 0
    else:
        fdr = len(fp)/(len(fp) + len(tp))

    tpr, fpr = len(tp)/len(H1_edges), len(fp)/len(H0_edges)
    return fpr, tpr, fdr, precision, recall


def run_lasso(data_1, data_2, real_diff_nodes, ols_delta_hat, c=30):
    s_cov_1 = np.cov(data_1.T)
    s_cov_2 = np.cov(data_2.T)
    evals = []
    e_values = []

    r_t0 = time.perf_counter()
    temp = SPDtrace(s_cov_2, s_cov_1, sparsityLevel=c, verbose=False)
    active_sets = []
    for i in range(len(temp["solution_path"])):
        for j in range(len(temp["solution_path"][i]["active_set"])):
            active_sets.append((int(temp["solution_path"][i]["active_set"][j]),
            int(temp["solution_path"][i]["active_set_signs"][j])))
        estimated_graph = construct_graph(s_cov_1.shape[0], active_sets)
        mirror_stats = generate_mirror_statistics_sgn(ols_delta_hat, estimated_graph)
        eval = np.max(get_e_value_nikahd(mirror_stats, 1))
        evals.append(eval)
        e_values.append(get_e_value_nikahd(mirror_stats, 1))
    r_sec = time.perf_counter() - r_t0
    print(f'SPDtrace R algorithm is finished with time: {r_sec}')
    return e_values, evals


def get_e_value_nikahd(mirror_statistics, cutoff_value):
    sub_mat = list(mirror_statistics[np.triu_indices(d,k=0)])
    t1 = len(sub_mat)
    t2 = (mirror_statistics >= cutoff_value).astype(int)
    t3 = 1 + len([i for i in sub_mat if i <= -cutoff_value])
    return t1 * t2 / t3


def get_metrics_values(metric_hist):
    precision_hist, recall_hist, fpr_hist, tpr_hist, fdr_hist, q_hist, time_hist = [], [], [], [], [], [], []
    for metric_list in metric_hist:
        if not isinstance(metric_list, list):
            fpr_list = metric_list["FPR"].tolist()
            tpr_list = metric_list["TPR"].tolist()
            recall_list = metric_list["NRecall"].tolist()
            precision_list = metric_list["NPrecision"].tolist()
            fdr_hist = metric_list["FDR"].tolist()
            q_hist = []
            time_hist = []
        else:    
            if len(metric_list) == 0:
                continue
            metric_list = [t for t in metric_list if len(t) >= 6]
            precision_list = [t[3] for t in metric_list]
            recall_list = [t[4] for t in metric_list]
            fpr_list = [t[0] for t in metric_list]
            tpr_list = [t[1] for t in metric_list]
            fdr_hist.append([t[2] for t in metric_list])
            precision_hist.append(precision_list)
            recall_hist.append(recall_list)
            fpr_hist.append(fpr_list)
            tpr_hist.append(tpr_list)
            q_hist.append([t[5] for t in metric_list])
            time_hist.append([t[6] for t in metric_list])
        
    tpr_hist = np.array(list(filter(lambda x: x is not None, tpr_hist)))
    fpr_hist = np.array(list(filter(lambda x: x is not None, fpr_hist)))
    precision_hist = np.array(list(filter(lambda x: x is not None, precision_hist)))
    recall_hist = np.array(list(filter(lambda x: x is not None, recall_hist)))
    fdr_hist = np.array(list(filter(lambda x: x is not None, fdr_hist)))
    q_hist = np.array(list(filter(lambda x: x is not None, q_hist)))
    time_hist = np.array(list(filter(lambda x: x is not None, time_hist)))

    tpr_mean = column_mean_ignore_none(tpr_hist)
    fpr_mean = column_mean_ignore_none(fpr_hist)
    precision_mean = column_mean_ignore_none(precision_hist)
    recall_mean = column_mean_ignore_none(recall_hist)
    fdr_mean = column_mean_ignore_none(fdr_hist)
    q_mean = column_mean_ignore_none(q_hist)
    time_mean = column_mean_ignore_none(time_hist)
    time_sd = column_sd_ignore_none(time_hist)
    return tpr_mean, fpr_mean, precision_mean, recall_mean, fdr_mean, q_mean, time_mean, time_sd


def column_mean_ignore_none(array):
    means = []
    # Transpose the array to iterate over columns (or use axis=0 in masked operations)
    for col in array.T:
        # Mask the None values in each column
        masked_col = np.ma.masked_equal(col, None)

        # Check if all values in the column are None
        if masked_col.count() == 0:
            means.append(np.nan)  # Append NaN if no valid values in the column
        else:
            means.append(masked_col.mean())  # Calculate mean for columns with valid values

    return np.array(means)

def column_sd_ignore_none(array):
    sds = []
    # Transpose the array to iterate over columns (or use axis=0 in masked operations)
    for col in array.T:
        # Mask the None values in each column
        masked_col = np.ma.masked_equal(col, None)

        # Check if all values in the column are None
        if masked_col.count() == 0:
            sds.append(np.nan)  # Append NaN if no valid values in the column
        else:
            sds.append(masked_col.std(ddof=1))  # Calculate standard deviation for columns with valid values

    return np.array(sds)



def extract_metrics(metric_list):
    precision_list = []
    recall_list = []
    fdr_list = []
    qs_list = []
    time_list_mean = []
    time_list_sd = []
    for i in range(len(metric_list)):
        metric_list_hist = metric_list[i]
        tpr_hist_mean, fpr_hist_mean, precision_hist_mean, recall_hist_mean, fdr_hist_mean, qs_hist_mean, time_hist_mean, time_hist_sd = get_metrics_values(metric_list_hist)
        precision_list.append(precision_hist_mean)
        recall_list.append(recall_hist_mean)
        fdr_list.append(fdr_hist_mean)
        qs_list.append(qs_hist_mean)
        time_list_mean.append(time_hist_mean)
        time_list_sd.append(time_hist_sd)
    metrics = [precision_list, recall_list, fdr_list, qs_list, time_list_mean, time_list_sd]
    return metrics


def save_metric_files(metric_1):
    with open("Data/our_method_metric_list", "wb") as fp:
       pickle.dump(metric_1, fp)


# Sign based Mirror Statistic with different nodes
d_list = [100, 200]
n = 100
s = 10
c = 15
rep = 50
dimension_metric_list = []
dimension_metric_list_DNetFinder = []
dimension_metric_list_DiffNetFDR_Xia2015 = []
dimension_metric_list_DiffNetFDR_Liu2017 = []

for d in d_list:
    sgn_metric_list_hist = []
    sgn_metric_list_hist_DNetFinder = []
    sgn_metric_list_hist_DiffNetFDR_Xia2015 = []
    sgn_metric_list_hist_DiffNetFDR_Liu2017 = []
    qs_list = [i/20 for i in range(2, 21)]

    for i in range(rep):
        print(f'dimension {d} run {i} ------------------------------------------------------------')
        # data generation
        dataset_1, precision_1, cov_1, dataset_2, precision_2, cov_2 = generate_reference_models(
            number_of_nodes=d, number_of_samples=n*d, number_of_changes=s, type="Full", mult=1)
        delta_star = precision_1 - precision_2
        real_H0, real_H1 = get_H0_H1(precision_1, precision_2)
        real_diff_nodes = set(get_diff_nodes(delta_star))


        # initial dataset split
        our_method_t0 = time.perf_counter()
        lasso_d1, lasso_d2, ols_d1, ols_d2 = split_dataset_samples(dataset_1, dataset_2, split_ratio=0.5)

        # get ols delta hat
        ols_delta_hat = get_delta_hat(ols_d1, ols_d2)

        # get lasso delta hat
        e_values, evals = run_lasso(lasso_d1, lasso_d2, real_diff_nodes, ols_delta_hat, c=c*s)
        our_method_sec = time.perf_counter() - our_method_t0
        p = len(list(ols_delta_hat[np.triu_indices(d,k=0)]))
        l = len(evals)
        print('Our algorithm is finised with time: ', our_method_sec)

        metric_list = []

        try:
            diffnetfdr_Xia2015_t0 = time.perf_counter()
            metric_list_DiffNetFDR_Xia2015 = DiffNetFDR_Xia2015(dataset_1, dataset_2, qs_list, delta_star)
            diffnetfdr_Xia2015_sec = time.perf_counter() - diffnetfdr_Xia2015_t0
            print(f'diffnetfdr_Xia2015 algorithm is finished with time: {diffnetfdr_Xia2015_sec}')
            metric_list_DiffNetFDR_Xia2015["diffnetfinder_sec"] = diffnetfdr_Xia2015_sec
            metric_list_DiffNetFDR_Xia2015['run'] = i
            metric_list_DiffNetFDR_Xia2015['dim'] = d
            sgn_metric_list_hist_DiffNetFDR_Xia2015.append(metric_list_DiffNetFDR_Xia2015)

        except Exception as e:
            diffnetfinder_sec = float("nan")
            metric_list_DiffNetFDR_Xia2015 = pd.DataFrame(
                {
                    "alpha": qs_list,
                    "FPR": [np.nan] * len(qs_list),
                    "TPR": [np.nan] * len(qs_list),
                    "FDR": [np.nan] * len(qs_list),
                    "NPrecision": [np.nan] * len(qs_list),
                    "NRecall": [np.nan] * len(qs_list),
                    "diffnetfinder_sec": [diffnetfinder_sec] * len(qs_list),
                    "run": [i] * len(qs_list),
                    "dim": [d] * len(qs_list),
                }
            )
            sgn_metric_list_hist_DiffNetFDR_Xia2015.append(metric_list_DiffNetFDR_Xia2015)
            print(i, 'sgn_metric_list_hist_DiffNetFDR_Xia2015')

        try:
            diffnetfdr_Liu2017_t0 = time.perf_counter()
            metric_list_DiffNetFDR_Liu2017 = DiffNetFDR_Liu2017(dataset_1, dataset_2, qs_list, delta_star)
            diffnetfdr_Liu2017_sec = time.perf_counter() - diffnetfdr_Liu2017_t0
            print(f'diffnetfdr_Liu2017 algorithm is finished with time: {diffnetfdr_Liu2017_sec}')
            metric_list_DiffNetFDR_Liu2017["diffnetfinder_sec"] = diffnetfdr_Liu2017_sec
            metric_list_DiffNetFDR_Liu2017['run'] = i
            metric_list_DiffNetFDR_Liu2017['dim'] = d
            sgn_metric_list_hist_DiffNetFDR_Liu2017.append(metric_list_DiffNetFDR_Liu2017)

        except Exception as e:
            diffnetfinder_sec = float("nan")
            metric_list_DiffNetFDR_Liu2017 = pd.DataFrame(
                {
                    "alpha": qs_list,
                    "FPR": [np.nan] * len(qs_list),
                    "TPR": [np.nan] * len(qs_list),
                    "FDR": [np.nan] * len(qs_list),
                    "NPrecision": [np.nan] * len(qs_list),
                    "NRecall": [np.nan] * len(qs_list),
                    "diffnetfinder_sec": [diffnetfinder_sec] * len(qs_list),
                    "run": [i] * len(qs_list),
                    "dim": [d] * len(qs_list),
                }
            )
            sgn_metric_list_hist_DiffNetFDR_Liu2017.append(metric_list_DiffNetFDR_Liu2017)
            print(i, 'sgn_metric_list_hist_DiffNetFDR_Liu2017')

        # try:
        #     dnetfinder_t0 = time.perf_counter()
        #     metric_list_DNetFinder = DNetFinder_Liu2017(dataset_1, dataset_2, qs_list, delta_star)
        #     dnetfinder_sec = time.perf_counter() - dnetfinder_t0
        #     print(f'DNetFinder algorithm is finished with time: {dnetfinder_sec}')
        #     metric_list_DNetFinder["dnetfinder_sec"] = dnetfinder_sec
        #     metric_list_DNetFinder['run'] = i
        #     metric_list_DNetFinder['dim'] = d
        #     sgn_metric_list_hist_DNetFinder.append(metric_list_DNetFinder)

        # except Exception as e:
        #     dnetfinder_sec = float("nan")
        #     metric_list_DNetFinder = pd.DataFrame(
        #         {
        #             "alpha": qs_list,
        #             "FPR": [np.nan] * len(qs_list),
        #             "TPR": [np.nan] * len(qs_list),
        #             "FDR": [np.nan] * len(qs_list),
        #             "NPrecision": [np.nan] * len(qs_list),
        #             "NRecall": [np.nan] * len(qs_list),
        #             "dnetfinder_sec": [dnetfinder_sec] * len(qs_list),
        #             "run": [i] * len(qs_list),
        #             "dim": [d] * len(qs_list),
        #         }
        #     )
        #     sgn_metric_list_hist_DNetFinder.append(metric_list_DNetFinder)
        #     print(i, 'sgn_metric_list_hist_DNetFinder')

        for j1 in range(len(qs_list)):
            qs = qs_list[j1]
            k_hat = 0
            max_e = 0
            evals[0] = np.nan
            for k in range(1,l):
                e = np.count_nonzero(list(e_values[k][np.triu_indices(d,k=0)]))
                if e > 0 and evals[k] != 0 and evals[k] >= p / (qs * e) and max_e < e:
                    k_hat = k
                    max_e = e
            eval = np.max(evals[k_hat])
            if max_e == 0:
              eval = np.nan
            fpr, tpr, fdr, precision, recall = get_fpr_tpr(real_H0, real_H1, e_values[k_hat], evals[k_hat])
            metric_list.append((fpr, tpr, fdr, precision, recall, qs, our_method_sec))
        sgn_metric_list_hist.append(metric_list)

    dimension_metric_list.append(sgn_metric_list_hist)
    
    d_sgn_metric_list_hist_DiffNetFDR = pd.concat(sgn_metric_list_hist_DiffNetFDR_Xia2015, ignore_index=True)
    dimension_metric_list_DiffNetFDR_Xia2015.append(sgn_metric_list_hist_DiffNetFDR_Xia2015)
    d_sgn_metric_list_hist_DiffNetFDR.to_csv(f'Data/{d}_sgn_metric_list_hist_DiffNetFDR_Xia2015.csv', index=False)
    
    d_sgn_metric_list_hist_DiffNetFDR_Liu2017 = pd.concat(sgn_metric_list_hist_DiffNetFDR_Liu2017, ignore_index=True)
    dimension_metric_list_DiffNetFDR_Liu2017.append(sgn_metric_list_hist_DiffNetFDR_Liu2017)
    d_sgn_metric_list_hist_DiffNetFDR_Liu2017.to_csv(f'Data/{d}_sgn_metric_list_hist_DiffNetFDR_Liu2017.csv', index=False)

our_metrics = extract_metrics(dimension_metric_list)

save_metric_files(our_metrics)
print('Our Algorithm is finished.')
