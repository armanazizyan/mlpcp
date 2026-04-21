# library(RSNNS)
# library(foreach)
# library(doSNOW)

#' Two-sided moving average smoothing
#' @param x a numeric vector
#' @param n number of neighbors to smooth over in total
#' @param circular if the ends should wrap or remain NA
#' @return smoothed numeric vector
#' @importFrom stats filter
#' @export
ma <- function(x, n = 5, circular = T){
  stats::filter(x, rep(1 / n, n), sides = 2, circular = circular)
}


#' Fitting MLP for CP detection algorithm
#' @description
#' Fits MLPs with defined structure to rolling window of defined size,
#' Then fits another MLP with double the rolling window,
#' capturing two windows' uncertainty
#'
#' @param vec original data
#' @param act1 activation function of smaller MLP
#' @param act2 activation function of larger MLP
#' @param lr1,lr2 learning rates of smaller and larger MLP
#' @param hl1,hl2 number of neurons in the hidden layers of smaller and larger MLP
#' @param ep1,ep2 number of training epochs for smaller and larger MLP
#' @return a list of two elements containing the fitted values for the smaller
#' and larger MLP models, original input and index.
#' @import RSNNS foreach doSNOW parallel
#' @export
fit_mlp <- function(vec, w=100,
                 act1="Act_TanH_Xdiv2",
                 act2="Act_TanH_Xdiv2",
                 lr1=0.01, lr2=0.01, hl1=16,
                 hl2=32, ep1=1000, ep2=2000){
  n.val <- length(vec)
  # starts <- 1:(n.val-w+1)
  x <- 1:n.val
  y <- vec

  pb <- txtProgressBar(min = 0, max = n.val-w, style = 3)

  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  num_cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(num_cores)
  doSNOW::registerDoSNOW(cl)
  start_time <- Sys.time()

  res.list <- foreach(idx = 1:(n.val-w+1), .packages = c("RSNNS","mlpcp"),
                      .options.snow = opts) %dopar% {

    i <- idx

    sub_x <- as.matrix(as.vector(scale(x[i:(i + w-1)])))
    sub_y <- as.matrix(y[i:(i + w-1)])

    net <- mlp(
      sub_x, sub_y,
      size = hl1,
      learnFuncParams = lr1,
      maxit = ep1,
      linOut = TRUE,
      hiddenActFunc = act1 #"Act_TanH_Xdiv2" #

    )

    y_hat <- net$fitted.values
    cbind(sub_x, y_hat, i:(i + w-1))
  }


  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat("\n")
  print(time_taken)


  #stopCluster(cl)
  close(pb)


  # window_step <- 1
  # starts <- 1:(n.val-2*w+1)

  #num_cores <- parallel::detectCores() - 1
  #cl <- makeCluster(num_cores)
  #registerdoSNOW(cl)

  pb <- txtProgressBar(min = 0, max = n.val-2*w, style = 3)

  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  start_time <- Sys.time()

  res.list.dbl <- foreach(idx = 1:(n.val-2*w+1), .packages = c("RSNNS","mlpcp"),
                          .options.snow = opts) %dopar% {

    i <- idx


    sub_x <- as.matrix(as.vector(scale(x[i:(i + 2*w-1)])))
    sub_y <- as.matrix(y[i:(i + 2*w-1)])


    net <- mlp(
      sub_x, sub_y,
      size = hl2,
      learnFuncParams = lr2,
      maxit = ep2,
      linOut = TRUE,
      hiddenActFunc = act2#"Act_TanH_Xdiv2"#
    )

    y_hat <- net$fitted.values

    cbind(sub_x, y_hat, i:(i + 2*w-1))
  }

  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat("\n")
  print(time_taken)


  stopCluster(cl)

  close(pb)

  return(list(res.list, res.list.dbl))
}



#' Calculate Detector Statistic
#' @details
#' Calculates detector statistic using a linear combination of Ratio and DIfference
#' between RSS of full MLP and half window MLP
#'
#' @param y original data
#' @param fit_mlp_res output of fit_mlp function
#' @param w window size
#' @param a a numeric value between 0 and 1, 0 indicating only ratio, 1 indicating only difference
#' @param b correction to discontinuity for ratio
#' @param scale_01 should the detector be scaled to a 0 to 1 scale.
#' @returns a vector of detector statistic
#' @export
calc_detector <- function(y, fit_mlp_res, w=100, a=1, b=0, scale_01=T){
  res.list <- fit_mlp_res[[1]]
  res.list.dbl <- fit_mlp_res[[2]]
  n.val <- length(y)

  diff1 <- c()
  diff2 <- c()


  for(i in 1:(n.val-2*w)){
    rss1 <- sum((y[i:(i+w-1)]-res.list[[i]][,2])^2)
    rss2 <- sum((y[(i+w):(i+2*w-1)]-res.list[[i+w]][,2])^2)
    rss.tot <- sum((y[i:(i+2*w-1)]-res.list.dbl[[i]][,2])^2)
    diff1 <- c(diff1,((rss.tot+b)/(rss1+rss2+b)) )
  }

  for(i in 1:(n.val-2*w)){
    rss1 <- sum((y[i:(i+w-1)]-res.list[[i]][,2])^2)
    rss2 <- sum((y[(i+w):(i+2*w-1)]-res.list[[i+w]][,2])^2)
    rss.tot <- sum((y[i:(i+2*w-1)]-res.list.dbl[[i]][,2])^2)
    diff2 <- c(diff2,( (rss.tot-rss1-rss2)) )
  }

  if(scale_01){
    diff2 <- (diff2-min(diff2))/(max(diff2-min(diff2)))
    #
    diff1 <- (diff1-min(diff1))/(max(diff1-min(diff1)))

    diff <- (1-a)*diff1+a*diff2

    diff <- (diff-min(diff))/(max(diff-min(diff)))
  }
  else{
    diff <- (1-a)*diff1+a*diff2
  }

  diff
}

#' Best Split Free Assignment
#' @details
#' to be added
#'
#'
#' @export
#'
best_split_free <- function(labels) {
  # Validate
  if (length(labels) < 2) stop("Need at least 2 labels to split.")
  if (any(is.na(labels))) stop("Labels contain NA; please remove or impute.")
  # Normalize labels to 1/2 as integers
  u <- sort(unique(labels))
  if (!all(u %in% c(1, 2))) stop("Labels must be in {1, 2}.")

  n <- length(labels)
  # cumulative counts
  cum1 <- cumsum(labels == 1)
  cum2 <- cumsum(labels == 2)
  total1 <- cum1[n]
  total2 <- cum2[n]

  # For cut at i (split between i and i+1), i in 1..(n-1)
  # Left correct = max(#1_left, #2_left)
  # Right correct = max(#1_right, #2_right)
  left_correct  <- pmax(cum1, cum2)
  right_correct <- pmax(total1 - cum1, total2 - cum2)
  total_correct <- left_correct + right_correct

  # Only consider cuts 1..(n-1)
  total_correct <- total_correct[1:(n - 1)]

  # If multiple maxima, pick the earliest cut (or change to which.max(last) for latest)
  best_i <- which.max(total_correct)

  # Return a small list with details
  list(
    index = best_i,                             # split is between index and index+1
    total_correct = total_correct[best_i],      # max number correctly separated
    left_majority  = if (cum1[best_i] >= cum2[best_i]) 1 else 2,
    right_majority = if ((total1 - cum1[best_i]) >= (total2 - cum2[best_i])) 1 else 2,
    accuracy_if_assigning_majorities = total_correct[best_i] / n
  )
}

#' Best Spike Select
#' @details
#' to be added
#'
#'
#' @export
#'
select_best_spike <- function(s, right_tail_cutoff = 0.95,
                              left_tail_cutoff = 0.6) {
  # s: smoothed spacing vector
  # p: ECDF probabilities (same length as s)
  # tail_cutoff: ignore spikes with p > tail_cutoff

  p <- (1:length(s)) / (length(s) + 1)


  # Find peaks using pracma::findpeaks
  peaks <- pracma::findpeaks(s, sortstr = FALSE)
  if (is.null(peaks)) return(NA)  # no peaks found

  # Extract peak components
  peak_heights <- peaks[, 1]
  peak_idx     <- peaks[, 2]
  left_base    <- peaks[, 3]
  right_base   <- peaks[, 4]

  # Compute prominence
  baseline <- pmax(s[left_base], s[right_base])
  prominence <- peak_heights - baseline

  # Exclude tail region (avoid trivial spike near p=1)
  valid <- (p[peak_idx] < right_tail_cutoff) & (p[peak_idx] > left_tail_cutoff)
  if (!any(valid)) return(NA)

  # Select the spike with largest prominence
  best_idx <- peak_idx[valid][which.max(prominence[valid])]
  best_p   <- p[best_idx]

  return(best_p)
}


#' Signal Decomposition
#' @details
#' to be added
#'
#'
#' @export
#'
decompose_signal <- function(
    y,
    diff,
    w = 200,
    ma_window = 100,
    right_tail_cutoff = 0.95,
    left_tail_cutoff = 0.6,
    threshold = "auto"
) {

  n <- length(y)

  ## 1. Local maxima of raw detector
  pre.sm.local <- pracma::findpeaks(diff, minpeakdistance = w)

  ## 2. Smooth detector
  ma.diff <- as.numeric(ma(diff, ma_window))

  ## 3. Local maxima of smoothed detector
  l.max <- pracma::findpeaks(ma.diff, minpeakdistance = w)

  ## 4. ECDF significance
  f_diff <- ecdf(ma.diff)
  ecdf_vals <- f_diff(l.max[,1])
  l.max <- cbind(l.max, ecdf_vals)  # col 5 = ECDF

  ## 5. ECDF spacing curve
  x.diff <- sort(ma.diff)
  dx <- diff(x.diff)
  s <- as.numeric(ma(dx, w, circular = FALSE))
  s[is.na(s)] <- 0

  ## 6. Threshold selection
  if (is.character(threshold) && threshold == "auto") {
    thr <- select_best_spike(
      s,
      right_tail_cutoff = right_tail_cutoff,
      left_tail_cutoff  = left_tail_cutoff
    )
  } else if (is.numeric(threshold)) {
    thr <- threshold
  } else {
    stop("threshold must be 'auto' or numeric.")
  }

  ## 7. Changepoints
  cps_raw <- l.max[l.max[,5] >= thr, , drop = FALSE]
  cps <- cps_raw[,2] + w

  ## 8. Local correction around each changepoint
  marg <- w
  cor.cps <- c()
  shift.vals <- c()
  wc.p <- c()

  for (cp in sort(cps)) {

    lo <- max(1, cp - marg)
    hi <- min(n, cp + marg)

    segment <- y[lo:hi]

    b <- kmeans(segment, centers = 2)
    bsf <- best_split_free(b$cluster)

    corrected_cp <- cp + (marg - bsf$index)
    corrected_cp <- min(max(corrected_cp, 1), n)
    cor.cps <- c(cor.cps, corrected_cp)

    # safe windows
    w1_start <- max(1, cp - bsf$index)
    w1_end   <- min(n, cp + marg - bsf$index)

    w2_start <- max(1, (cp + marg - bsf$index))
    w2_end   <- min(n, (cp + 2*marg - bsf$index))

    shift.vals <- c(
      shift.vals,
      median(y[w1_start:w1_end]) - median(y[w2_start:w2_end])
    )

    wc.p <- c(
      wc.p,
      wilcox.test(
        y[w1_start:w1_end],
        y[w2_start:w2_end],
        exact = FALSE
      )$p.value
    )
  }

  ## 9. Raw piecewise correction vector (cor.vec)
  cps <- cor.cps
  starts <- cps + 1
  ends <- c(cps[-1], n)

  cor.vec <- numeric(n)
  cumulative.shift.vals <- cumsum(shift.vals)

  for (i in seq_along(starts)) {
    cor.vec[starts[i]:ends[i]] <- cumulative.shift.vals[i]
  }

  ## 10. Apply correction
  new.y <- y + cor.vec

  ## 11. MLP fit on corrected signal
  x <- 1:n
  sub_x <- as.matrix(as.vector(scale(x)))
  sub_y <- as.matrix(new.y)

  net <- RSNNS::mlp(
    sub_x, sub_y,
    size = 64,
    learnFuncParams = c(0.01, 0.001),
    maxit = 2000,
    linOut = TRUE,
    hiddenActFunc = "Act_Logistic"
  )

  new.y_hat <- net$fitted.values
  new.resid <- new.y - new.y_hat

  ## 12. Step means of residuals (corrected piecewise constant)
  idx_start <- c(1, cps + 1)
  idx_end <- c(cps, length(y))

  step_mean <- numeric(length(y))
  for (i in seq_along(idx_start)) {
    step_mean[idx_start[i]:idx_end[i]] <-
      mean((y - new.y_hat)[idx_start[i]:idx_end[i]])
  }

  min.step <- step_mean[1]
  step_mean <- step_mean - min.step
  new.y_hat <- new.y_hat + min.step

  ## 13. Return
  list(
    corrected_signal   = new.y,
    smooth_curve       = new.y_hat,
    piecewise_constant = step_mean,   # <- corrected piecewise constant
    raw_correction     = cor.vec,     # original cor.vec kept for reference
    residual           = new.resid,
    changepoints       = cps,
    changepoint_ecdf   = cps_raw[,5],
    shift_values       = shift.vals,
    p_values           = wc.p,
    smoothed_detector  = ma.diff,
    ecdf_spacing       = s,
    threshold_used     = thr,
    local_maxima       = l.max,
    mlp_model          = net
  )
}


#' Confidence Interval
#' @details
#' to be added
#'
#'
#' @export
#'
changepoint_CI <- function(
    decomp,
    alpha = 0.05,
    B = 5000,
    K = 500,
    regime = c("auto", "vanishing", "nonvanishing")
) {
  regime <- match.arg(regime)

  cps        <- decomp$changepoints
  xi.hat.vec <- decomp$shift_values
  resid      <- decomp$residual
  n          <- length(resid)

  if (length(cps) == 0L) {
    stop("No changepoints found in 'decomp'.")
  }

  # noise estimate
  sigma.hat <- sd(resid)

  # helper: simulate vanishing-jump limit (Brownian with drift)
  simulate_vanishing <- function(B, K, sigma) {
    Zmax <- numeric(B)
    zeta <- -K:K
    m    <- length(zeta)
    for (b in seq_len(B)) {
      W <- cumsum(rnorm(m, 0, 1))
      Z <- 2 * sigma * W - abs(zeta)
      Zmax[b] <- zeta[which.max(Z)]
    }
    Zmax
  }

  # helper: simulate non-vanishing-jump limit (random walk with drift)
  simulate_nonvanishing <- function(B, K, xi, sigma) {
    Zmax <- numeric(B)
    zeta <- -K:K
    for (b in seq_len(B)) {
      zpos <- rnorm(K, mean = -xi^2, sd = sqrt(4 * xi^2 * sigma^2))
      zneg <- rnorm(K, mean = -xi^2, sd = sqrt(4 * xi^2 * sigma^2))

      Cpos <- cumsum(zpos)
      Cneg <- cumsum(zneg)

      C <- c(rev(-Cneg), 0, Cpos)
      Zmax[b] <- zeta[which.max(C)]
    }
    Zmax
  }

  # precompute Brownian limit once (parameter-free after scaling by xi^2)
  Z_brown <- simulate_vanishing(B = B, K = K, sigma = sigma.hat)
  q_brown <- quantile(Z_brown, c(alpha / 2, 1 - alpha / 2))

  # container
  CI_mat <- matrix(NA_real_, nrow = length(cps), ncol = 2)
  colnames(CI_mat) <- c("lower", "upper")

  for (j in seq_along(cps)) {
    tau.hat <- cps[j]
    xi.hat  <- xi.hat.vec[j]

    # decide regime
    reg_j <- regime
    if (regime == "auto") {
      reg_j <- if (xi.hat < sigma.hat) "vanishing" else "nonvanishing"
    }

    if (reg_j == "vanishing") {
      # xi^2 (tau_hat - tau_0) -> Z
      CI <- tau.hat + as.numeric(q_brown) / (xi.hat^2)
    } else {
      # simulate random-walk limit for this xi
      Z_rw <- simulate_nonvanishing(B = B, K = K, xi = xi.hat, sigma = sigma.hat)
      q_rw <- quantile(Z_rw, c(alpha / 2, 1 - alpha / 2))
      CI   <- tau.hat + as.numeric(q_rw)
    }

    # clip and round
    CI <- round(pmax(1, pmin(n, CI)))
    CI_mat[j, ] <- CI
  }

  data.frame(
    changepoint = cps,
    shift_value = xi.hat.vec,
    CI_lower    = CI_mat[, 1],
    CI_upper    = CI_mat[, 2]
  )
}

