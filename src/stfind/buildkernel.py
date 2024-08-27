# https://medium.com/analytics-vidhya/how-to-create-a-python-library-7d5aea80cc3f

import numpy as np
from sklearn.preprocessing import scale
from statsmodels.nonparametric.bandwidths import select_bandwidth
from statsmodels.sandbox.nonparametric import kernels

from .util import isPD, nearestPD

# from statsmodels.stats.correlation_tools import cov_nearest


def cal_bandwidth(expr, bw_method="silverman"):

    # https://www.statsmodels.org/dev/_modules/statsmodels/nonparametric/bandwidths.html#select_bandwidth
    bw = select_bandwidth(expr, bw=bw_method, kernel=kernels.Gaussian())

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html#scipy.stats.gaussian_kde
    # import scipy.stats as stats
    # from scipy.stats import gaussian_kde
    # bw = stats.gaussian_kde(expr, bw_method = bw_method).factor

    return bw


def bandwidth_select(expr, method="silverman"):

    tot_bw = []
    for i in range(expr.shape[0]):
        try:
            res_bw = cal_bandwidth(expr[i], bw_method=method)
            tot_bw.append(res_bw)
        except Exception as e:
            print(f"Gene {i} : {str(e)}")

    median_bw = np.median([bw for bw in tot_bw if not np.isnan(bw)])

    return median_bw


def cal_kernel(coord, bandwidth):

    pairwise_sq_dists = np.sum(
        (coord[:, np.newaxis, :] - coord[np.newaxis, :, :]) ** 2, axis=-1
    )
    kernelmat = np.exp(-1 * pairwise_sq_dists / bandwidth)

    return kernelmat


def spatialkernel(expr, coord, bandwidthtype="silverman", userbandwidth=None):

    # Standardization
    # print("# Scale the expression of each gene.")
    expr = scale(expr, axis=1)

    # Calculate bandwidth
    if userbandwidth is None:
        bandwidth = bandwidth_select(expr, method=bandwidthtype)
        print(f"# The bandwidth is {round(bandwidth, 4)} using {bandwidthtype}.")
    else:
        bandwidth = userbandwidth
        print(f"# The bandwidth is {round(bandwidth, 4)} (set by user).")

    # Calculate the kernel matrix using the bandwidth
    coord_normalized = scale(coord)
    kernelmat = cal_kernel(coord=coord_normalized, bandwidth=bandwidth)

    return kernelmat


def buildchol(skernel):

    if isPD(skernel) is False:
        # skernel_nearPD = cov_nearest(skernel) # too slow
        skernel_nearPD = nearestPD(skernel)
        chol_skernel = np.linalg.cholesky(skernel_nearPD)
    else:
        chol_skernel = np.linalg.cholesky(skernel)

    return chol_skernel
