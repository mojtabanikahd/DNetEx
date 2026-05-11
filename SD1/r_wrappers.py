import numpy as np
from typing import Any, Dict, Optional

from rpy2 import robjects
from rpy2.robjects import NULL, numpy2ri, pandas2ri
from rpy2.robjects.language import LangVector
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.vectors import ListVector


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

def DiffNetFDR_Xia2015(XA, XB, alphas, delta_star):
    pandas2ri.activate()
    numpy2ri.activate()
    # Convert Python lists or arrays to numpy arrays and ensure they are C-contiguous
    XA = np.ascontiguousarray(XA)
    XB = np.ascontiguousarray(XB)
    alphas = np.ascontiguousarray(alphas)
    delta_star = np.ascontiguousarray(delta_star)
  
    # Source your R code
    robjects.r.source('R_codes/Libraray.R')
    r_func = robjects.globalenv['DiffNetFDR_Xia2015']

    # Convert Python arrays to R matrices
    r_XA = robjects.r.matrix(XA, nrow=XA.shape[0], ncol=XA.shape[1])
    r_XB = robjects.r.matrix(XB, nrow=XB.shape[0], ncol=XB.shape[1])
    r_alphas = FloatVector(alphas)
    r_delta_star = robjects.r.matrix(delta_star, nrow=delta_star.shape[0], ncol=delta_star.shape[1])
  
    # Call the R function
    output = r_func(r_XA, r_XB, r_alphas, r_delta_star)
    output = pandas2ri.rpy2py(output)

    return output


def _names_as_str_list(nms_attr: Any) -> Optional[list[str]]:
    if nms_attr is None or nms_attr is NULL:
        return None
    names_arr = np.asarray(nms_attr, dtype=str)
    if names_arr.shape == ():
        return [str(names_arr.item())]
    return [str(x) for x in names_arr.ravel()]


def _r_lang_to_python_str(expr: LangVector) -> str:
    base = importr("base")
    return " ".join(str(line).strip() for line in base.deparse(expr))


def _r_to_python(obj: Any) -> Any:
    if obj is NULL:
        return None
    if isinstance(obj, np.ndarray):
        return np.asarray(obj)
    if isinstance(obj, LangVector):
        return _r_lang_to_python_str(expr=obj)
    if isinstance(obj, ListVector):
        names = _names_as_str_list(getattr(obj, "names", None))
        if names is None:
            return [_r_to_python(obj[i]) for i in range(len(obj))]
        return {n: _r_to_python(obj.rx2(n)) for n in names}
    try:
        return np.asarray(obj)
    except Exception:
        return obj


def SPDtrace(
    CovA: np.ndarray,
    CovB: np.ndarray,
    sparsityLevel: int,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run R ``SPDtrace::SPDtrace`` (lasso-penalized sparse differential network).

    Requires the SPDtrace package to be installed in the same R used by Python/rpy2.

    Args:
        CovA: Estimated covariance matrix for class A (2-D ``float``, square).
        CovB: Estimated covariance matrix for class B (same shape as ``CovA``).
        sparsityLevel: Matches R ``sparsityLevel``; optional cap on differential edges.
            Use ``None`` to pass R ``NULL`` (R then sets a default from ``choose(nrow(CovA), 2)``).
        verbose: Maps to R ``verbose`` logical.

    Returns:
        Mapping with the documented R ``$`` names:

        * ``solution_path``: nested Python-native structure derived from ``ListVector`` entries
          (dicts keyed by knot / active-set field names where R assigns names).
        * ``last_differential_network``: ``numpy.ndarray`` (adjacency).
        * ``lambda_sequence``: 1-D ``numpy.ndarray`` of knot ``lambda`` values.
        * ``call``: deparsed ``match.call()`` via ``base::deparse``; can be long when R inlines
          large matrix arguments.

        Example: from directory ``SD1``, two ``p=10`` identity-like PSD covariances and
        ``SPDtrace(CovA, CovB, sparsityLevel=None, verbose=False)`` matches the R usage
        pattern in ``?SPDtrace``.

    """
    pandas2ri.activate()
    numpy2ri.activate()

    CovA_arr = np.ascontiguousarray(CovA, dtype=float)
    CovB_arr = np.ascontiguousarray(CovB, dtype=float)
    if CovA_arr.ndim != 2 or CovB_arr.ndim != 2:
        raise ValueError("CovA and CovB must be 2-D arrays")
    if CovA_arr.shape != CovB_arr.shape:
        raise ValueError("CovA and CovB must have identical shapes")

    spd_pkg = importr("SPDtrace")
    nrow, ncol = CovA_arr.shape
    r_CovA = robjects.r.matrix(CovA_arr, nrow=nrow, ncol=ncol)
    r_CovB = robjects.r.matrix(CovB_arr, nrow=nrow, ncol=ncol)

    r_out = spd_pkg.SPDtrace(r_CovA, r_CovB, sparsityLevel, verbose=verbose)

    return {
        "solution_path": _r_to_python(r_out.rx2("solution_path")),
        "last_differential_network": np.asarray(
            r_out.rx2("last_differential_network"), dtype=float
        ),
        "lambda_sequence": np.asarray(
            r_out.rx2("lambda_sequence"), dtype=float
        ).ravel(),
        "call": _r_lang_to_python_str(expr=r_out.rx2("call")),
    }
