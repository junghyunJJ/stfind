
# bw.nrd0(expr[i, ])
# bw.ucv(expr[i, ], lower = 0.01, upper = 1)
# bw.SJ(expr[i, ], method = "dpi")

# kedd::h.ucv(expr[i, ])

bandwidth_select <- function(expr, method = "Silverman") {
  N <- dim(expr)[2]
  if (method == "SJ") {
    bw_SJ <- c()
    for (i in 1:dim(expr)[1]) { # nolint: seq_linter.
      tryCatch(
        {
          # print(i)
          bw_SJ[i] <- bw.SJ(expr[i, ], method = "dpi")
        },
        error = function(e) {
          cat("Gene", i, " :", conditionMessage(e), "\n")
        }
      )
    }

    beta <- median(na.omit(bw_SJ))
  } else if (method == "Silverman") {
    bw_Silverman <- c()
    for (i in 1:dim(expr)[1]) { # nolint: seq_linter.
      tryCatch(
        {
          bw_Silverman[i] <- bw.nrd0(expr[i, ])
        },
        error = function(e) {
          cat("Gene", i, " :", conditionMessage(e), "\n")
        }
      )
    }
    beta <- median(na.omit(bw_Silverman))
  } else if (method == "mlcv") {
    bw_mlcv <- c()
    for (i in 1:dim(expr)[1]) { # nolint: seq_linter.
      tryCatch(
        {
          bw_mlcv[i] <- kedd::h.mlcv(expr[i, ]) # # https://cran.r-project.org/web/packages/kedd/vignettes/kedd.pdf
        },
        error = function(e) {
          cat("Gene", i, " :", conditionMessage(e), "\n")
        }
      )
    }
    beta <- median(na.omit(bw_mlcv))
  }
}

cal_kernel <- function(coord, bandwidth) {
  pairwise_sq_dists <- as.matrix(dist(coord)^2)
  K <- exp(-1 * pairwise_sq_dists / bandwidth)
  return(K)
}


spatialkernel <- function(exp, coord, bandwidthtype = "Silverman", userbandwidth = NULL, ncores = 1) {

  # Please check the deim for exp and coord
  stopifnot()

  # 1. standadization
  cat(paste0("## Scale the expression of each gene. \n"))
  expr <- t(scale(t(exp)))

  # 2. cal bandwidth
  if (is.null(userbandwidth)) {
    bandwidth <- bandwidth_select(expr, method = bandwidthtype)
    cat(paste0("## The bandwidth is: ", round(bandwidth, 3), " using", bandwidthtype, " \n"))
  } else {
    bandwidth <- userbandwidth
    cat(paste0("## The bandwidth is: ", bandwidth, " (set by user)\n"))
  }

  # 3. Calculate the kernel matrix using the bandwidth
  coord_normalized <- scale(coord)
  kernelmat <- cal_kernel(coord = coord_normalized, bandwidth = bandwidth)

  return(kernelmat)
}

# exp <- data.table::fread("data/expr.csv")
# coord <- data.table::fread("data/coord.csv")

# # spatialkernel(exp, coord)[1:5, 1:5]
# spatialkernel(exp, coord, userbandwidth = 0.1423)[1:5, 1:5]
