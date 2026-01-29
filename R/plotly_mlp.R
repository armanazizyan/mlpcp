#library(plotly)


#' Interactive plot of fit_mlp output
#' @param y the original vector of observations
#' @param fit_mlp_res the resulting list of fit_mlp function
#' @return plotly object
#' @import plotly
#' @export
plot_mlp_slider_3lines <- function(
    y, fit_mlp_res,
    t.chp.ind=NA, step = 50, start = 75, w = 100
) {

  n <- length(y)

  res.list1 <- fit_mlp_res[[2]]
  res.list2 <- fit_mlp_res[[1]]
  res.list3 <- fit_mlp_res[[1]][w+1:length(y)]

  idx_list <- seq(start, n - 2*w, by = step)
  n_models <- length(idx_list)

  # --- CP segments ---
  if(!is.na(t.chp.ind)){
    cp_segments <- lapply(t.chp.ind, function(cp) {
      list(
        type = "line",
        x0 = cp, x1 = cp,
        y0 = min(y), y1 = max(y),
        line = list(color = "blue", dash = "dash")
      )
    })
  }else{
    cp_segments <- NA}


  # --- Initial traces (using res.list1, res.list2, res.list3) ---
  i0 <- idx_list[1]

  p <- plot_ly() %>%
    # trace 1: original y (static)
    add_lines(
      x = 1:n, y = y,
      name = "Original Signal",
      line = list(color = "rgba(0,0,0,0.3)")
    ) %>%
    # trace 2: model 1 (dynamic)
    add_lines(
      x = res.list1[[i0]][,3],
      y = res.list1[[i0]][,2],
      name = "Model 1",
      line = list(color = "red", width = 3)
    ) %>%
    # trace 3: model 2 (dynamic)
    add_lines(
      x = res.list2[[i0]][,3],
      y = res.list2[[i0]][,2],
      name = "Model 2",
      line = list(color = "green", width = 3)
    ) %>%
    # trace 4: model 3 (dynamic)
    add_lines(
      x = res.list3[[i0]][,3],
      y = res.list3[[i0]][,2],
      name = "Model 3",
      line = list(color = "purple", width = 3)
    ) %>%
    layout(
      title = "Multiple MLP Fits – Interactive Slider",
      xaxis = list(title = "Index"),
      yaxis = list(title = "Value"),
      shapes = cp_segments,

      sliders = list(
        list(
          active = 0,
          currentvalue = list(prefix = "Model index: "),
          steps = lapply(1:n_models, function(k) {
            list(
              label = as.character(idx_list[k]),
              method = "animate",
              args = list(
                list(paste0("frame", k)),
                list(mode = "immediate",
                     frame = list(duration = 0, redraw = FALSE),
                     transition = list(duration = 0))
              )
            )
          })
        )
      )
    )


  # --- Frames: update traces 2, 3, and 4 only ---
  frames <- lapply(1:n_models, function(k) {
    i <- idx_list[k]

    list(
      name = paste0("frame", k),
      data = list(
        NULL,  # trace 1 stays unchanged
        list(  # trace 2 (model 1)
          x = res.list1[[i]][,3],
          y = res.list1[[i]][,2],
          mode = "lines",
          line = list(color = "red", width = 3)
        ),
        list(  # trace 3 (model 2)
          x = res.list2[[i]][,3],
          y = res.list2[[i]][,2],
          mode = "lines",
          line = list(color = "green", width = 3)
        ),
        list(  # trace 4 (model 3)
          x = res.list3[[i]][,3],
          y = res.list3[[i]][,2],
          mode = "lines",
          line = list(color = "purple", width = 3)
        )
      )
    )
  })

  p$x$frames <- frames
  p

}
