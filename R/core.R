# library(RSNNS)
# library(foreach)
# library(doParallel)

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


#' Fits MLP with defined structure to rolling window of defined size,
#' Then fits another MLP with double the rolling window,
#' capturing two windows' uncertainty
#' @param vec original data
#' @param act1 activation function of smaller MLP
#' @param act2 activation function of larger MLP
#' @param lr1,lr2 learning rates of smaller and larger MLP
#' @param hl1,hl2 number of neurons in the hidden layers of smaller and larger MLP
#' @param ep1,ep2 number of training epochs for smaller and larger MLP
#' @return a list of two elements containing the fitted values for the smaller
#' and larger MLP models, original input and index.
#' @import RSNNS foreach doParallel parallel
#' @export
fit_mlp <- function(vec, w=100,
                 act1="Act_TanH_Xdiv2",
                 act2="Act_TanH_Xdiv2",
                 lr1=0.01, lr2=0.01, hl1=16,
                 hl2=32, ep1=1000, ep2=2000){
  n.val <- length(vec)
  starts <- 1:(n.val-w)
  x <- 1:n.val
  y <- vec

  num_cores <- parallel::detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  start_time <- Sys.time()

  res.list <- foreach(idx = seq_along(starts), .packages = "RSNNS") %dopar% {

    i <- starts[idx]

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
  print(time_taken)


  #stopCluster(cl)


  window_step <- 1
  starts <- 1:(n.val-w)

  num_cores <- parallel::detectCores() - 1
  #cl <- makeCluster(num_cores)
  #registerDoParallel(cl)


  start_time <- Sys.time()

  res.list.dbl <- foreach(idx = seq_along(starts), .packages = "RSNNS") %dopar% {

    i <- starts[idx]


    sub_x <- as.matrix(as.vector(scale(x[i:(i + w-1)])))
    sub_y <- as.matrix(y[i:(i + w-1)])


    net <- mlp(
      sub_x, sub_y,
      size = hl2,
      learnFuncParams = lr2,
      maxit = ep2,
      linOut = TRUE,
      hiddenActFunc = act2#"Act_TanH_Xdiv2"#
    )

    y_hat <- net$fitted.values

    cbind(sub_x, y_hat, i:(i + w-1))
  }

  end_time <- Sys.time()
  time_taken <- end_time - start_time
  print(time_taken)


  stopCluster(cl)

  return(list(res.list, res.list.dbl))
}

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

  diff2 <- (diff2-min(diff2))/(max(diff2)-min(diff2))

  diff1 <- (diff1-min(diff1))/(max(diff1)-min(diff1))

  diff <- (1-a)*diff1+a*diff2

  diff <- (diff-min(diff))/(max(diff)-min(diff))

  diff
}
