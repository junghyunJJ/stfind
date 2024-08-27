import os
import shutil

import numpy as np
import pandas as pd
from ctmm import log, op
from statsmodels.stats.multitest import fdrcorrection

from .buildkernel import buildchol, spatialkernel

# import multiprocessing as mp


def stfind(
    exprs,
    prop,
    coord=None,
    genes=None,
    bandwidthtype="silverman",
    tmpdir="./tmp",
    remove_tmpdir=True,
):
    log.logger.info("Start")

    # TODO: we need to check index names or colnames of 'exprs',' prop', and 'coord'

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
        os.makedirs(tmpdir)
    else:
        os.makedirs(tmpdir)

    n_genes, n_cells = exprs.shape
    celltypes = prop.columns

    # save celltype prop matrix
    P_f = f"{tmpdir}/prop.gz"
    prop.to_csv(P_f, sep="\t", index=False, header=False, compression="gzip")

    # save NULL nu
    nu_f = f"{tmpdir}/nu.gz"
    pd.DataFrame(np.repeat(0, n_cells)).to_csv(
        nu_f, sep="\t", index=False, header=False, compression="gzip"
    )

    # cal spatial kernel for 'stfind + spatial'
    if coord is not None:
        log.logger.info("Build a spatial kernel")
        skernel = spatialkernel(
            np.array(exprs), np.array(coord), bandwidthtype=bandwidthtype
        )
        chol_skernel = buildchol(skernel)

        chol_skernel_f = f"{tmpdir}/chol_skernel.gz"
        pd.DataFrame(chol_skernel).to_csv(
            chol_skernel_f, sep="\t", index=False, header=False, compression="gzip"
        )

    # filtering exprs if we provide the specific gene list
    if genes is not None:
        exprs = exprs[exprs.index.isin(genes)]
        n_genes, n_cells = exprs.shape

    log.logger.info(f"{n_cells} cells / {n_genes} genes/ {len(celltypes)} celltypes")

    save_res = list()

    # TODO!!! we need to update the for loop using multiprocessing?
    for i in range(len(exprs.index)):

        sel_gene = exprs.index[i]
        # TODO!!! we need to update the log output
        log.logger.info(f"{i+1}/{len(exprs.index)}: {sel_gene}")

        y_f = f"{tmpdir}/{sel_gene}_exp.gz"
        exprs.iloc[i].to_csv(
            y_f, sep="\t", index=False, header=False, compression="gzip"
        )

        if coord is None:
            log.logger.info("Run stfind")
            free, p_free = op.free_REML(
                y_f=y_f, P_f=P_f, nu_f=nu_f, method="BFGS", optim_by_R=True
            )
        elif coord is not None:
            log.logger.info("Run stfind with spatial kernel")
            free, p_free = op.free_REML(
                y_f=y_f,
                P_f=P_f,
                nu_f=nu_f,
                random_covars_d={"spatial": chol_skernel_f},
                method="BFGS",
                optim_by_R=True,
            )

        vc = pd.DataFrame(
            {
                "celltype": celltypes,
                "vc": np.diag(free["V"]),
                "p": p_free["Vi"],
                "fdr": fdrcorrection(p_free["Vi"])[1],
            }
        )
        res = {
            "vc": vc,
            "vs": free["randomeffect_vars"],
            "ve:": free["hom2"],
            "ll": free["l"],
            "r2": free["r2"],
        }
        save_res.append(res)

    # remove tmpdir
    if remove_tmpdir:
        shutil.rmtree(tmpdir)
    log.logger.info("End")
    return save_res
