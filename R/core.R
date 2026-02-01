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
#' @returns a vector of detector statistic
#' @export
calc_detector <- function(y, fit_mlp_res, w=100, a=1, b=0){
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

  diff2 <- (diff2-min(diff2))/(max(diff2-min(diff2)))
  #
  diff1 <- (diff1-min(diff1))/(max(diff1-min(diff1)))

  diff <- (1-a)*diff1+a*diff2

  diff <- (diff-min(diff))/(max(diff-min(diff)))

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
