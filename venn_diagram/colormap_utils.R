get_color_palette <- function(cmap, n = 200, circos = TRUE) {
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install the 'viridis' package.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Please install the 'RColorBrewer' package.")
  library(viridis)
  library(RColorBrewer)
  library(scico)

  reverse <- FALSE
  if (grepl("_r$", cmap)) {
    reverse <- TRUE
    cmap <- sub("_r$", "", cmap)
  }

  # 离散型调色板
  brewer_maps <- c(
    "PiYG","PRGn","BrBG","PuOr","RdGy","RdBu","RdYlBu","RdYlGn","Spectral",
    "YlGnBu","YlGn","YlOrBr","YlOrRd","GnBu","BuGn","BuPu","Purples","Greens",
    "Oranges","Reds","Greys","PuRd","RdPu","Blues",
    "Pastel1","Pastel2","Set1","Set2","Set3","Dark2","Accent",
    "Paired","OrRd","PuBu","PuBuGn"
  )
  brewer_discrete <- c(
  "Pastel1","Pastel2","Set1","Set2","Set3","Dark2","Accent",
  "Paired"
)
  brewer_continuous <- setdiff(brewer_maps, brewer_discrete)


  matplotlib_categorical <- list(
    tab10 = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
              "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"),
    tab20 = c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
              "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94",
              "#e377c2","#f7b6d2","#7f7f7f","#c7c7c7","#bcbd22","#dbdb8d",
              "#17becf","#9edae5"),
    tab20b= c("#393b79","#5254a3","#6b6ecf","#9c9ede",
           "#637939","#8ca252","#b5cf6b","#cedb9c",
           "#8c6d31","#bd9e39","#e7ba52","#e7cb94",
           "#843c39","#ad494a","#d6616b","#e7969c",
           "#7b4173","#a55194","#ce6dbd","#de9ed6"),
    tab20c= c("#3182bd","#6baed6","#9ecae1","#c6dbef",
           "#e6550d","#fd8d3c","#fdae6b","#fdd0a2",
           "#31a354","#74c476","#a1d99b","#c7e9c0",
           "#756bb1","#9e9ac8","#bcbddc","#dadaeb",
           "#636363","#969696","#bdbdbd","#d9d9d9")
  )

  viridis_maps <- c("viridis","plasma","inferno","magma","cividis")
  custom_maps <- c("coolwarm","bwr","seismic")
  scico_maps <- c("berlin","managua","vanimo")

  pal <- NULL

  # matplotlib 离散型
  if (cmap %in% names(matplotlib_categorical)) {
    base_colors <- matplotlib_categorical[[cmap]]
    if (circos) {
      pal <- rep(base_colors, length.out = n)
    } else {
      pal <- if (n <= length(base_colors)) base_colors[1:n] else colorRampPalette(base_colors)(n)
    }

  # RColorBrewer 离散型
  }else if (cmap %in% brewer_discrete) {
  max_colors <- RColorBrewer::brewer.pal.info[cmap, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(max_colors, cmap)
  if (circos) {
    pal <- rep(base_colors, length.out = n)
  } else {
    pal <- if (n <= max_colors) base_colors[1:n] else colorRampPalette(base_colors)(n)
  }

# RColorBrewer 连续型
} else if (cmap %in% brewer_continuous) {
  max_colors <- RColorBrewer::brewer.pal.info[cmap, "maxcolors"]
  base_colors <- RColorBrewer::brewer.pal(max_colors, cmap)
  pal <- colorRampPalette(base_colors)(n)

# viridis / custom / scico
}
            else if (cmap %in% viridis_maps) {
    pal <- switch(cmap,
                  viridis = viridis::viridis(n),
                  plasma  = viridis::plasma(n),
                  inferno = viridis::inferno(n),
                  magma   = viridis::magma(n),
                  cividis = viridis::cividis(n))

  # 自定义渐变
  } else if (cmap %in% custom_maps) {
    if (cmap == "coolwarm") pal <- colorRampPalette(c("#3B4CC0","#FFFFFF","#B40426"))(n)
    if (cmap == "bwr")      pal <- colorRampPalette(c("#0000FF","#FFFFFF","#FF0000"))(n)
    if (cmap == "seismic")  pal <- colorRampPalette(c("#67001F","#FFFFFF","#053061"))(n)

  # scico
  } else if (cmap %in% scico_maps) {
    pal <- scico::scico(n, palette = cmap)
  } else {
    stop("Unsupported colormap: ", cmap)
  }

  if (reverse) pal <- rev(pal)
  return(pal)
}
get_color_palette_old1 <- function(cmap, n = 200, circos = TRUE) {
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install 'viridis'.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Please install 'RColorBrewer'.")
  if (!requireNamespace("scico", quietly = TRUE)) stop("Please install 'scico'.")

  library(viridis)
  library(RColorBrewer)
  library(scico)

  reverse <- FALSE
  if (grepl("_r$", cmap)) {
    reverse <- TRUE
    cmap <- sub("_r$", "", cmap)
  }

  # 只保留 RColorBrewer 里的离散型 (qual)
  brewer_qual <- rownames(RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual", ])

  # matplotlib 离散型
  matplotlib_categorical <- list(
    tab10 = c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
              "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"),
    tab20 = c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
              "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94",
              "#e377c2","#f7b6d2","#7f7f7f","#c7c7c7","#bcbd22","#dbdb8d",
              "#17becf","#9edae5")
  )

  # 连续型调色板
  viridis_maps <- c("viridis","plasma","inferno","magma","cividis")
  custom_maps  <- c("coolwarm","bwr","seismic")
  scico_maps   <- c("berlin","managua","vanimo")

  pal <- NULL
  is_discrete <- FALSE

  ## 1) matplotlib 离散
  if (cmap %in% names(matplotlib_categorical)) {
    base_colors <- matplotlib_categorical[[cmap]]
    is_discrete <- TRUE
    if (circos) {
      pal <- rep(base_colors, length.out = n)
    } else {
      pal <- if (n <= length(base_colors)) base_colors[1:n] else colorRampPalette(base_colors)(n)
    }

  ## 2) RColorBrewer 离散型
  } else if (cmap %in% brewer_qual) {
    max_colors <- RColorBrewer::brewer.pal.info[cmap, "maxcolors"]
    base_colors <- RColorBrewer::brewer.pal(max_colors, cmap)
    is_discrete <- TRUE
    if (circos) {
      pal <- rep(base_colors, length.out = n)
    } else {
      if (n <= max_colors) {
        pal <- base_colors[1:n]
      } else {
        pal <- colorRampPalette(base_colors)(n)
      }
    }

  ## 3) 连续型 viridis
  } else if (cmap %in% viridis_maps) {
    pal <- switch(cmap,
                  viridis = viridis::viridis(n),
                  plasma  = viridis::plasma(n),
                  inferno = viridis::inferno(n),
                  magma   = viridis::magma(n),
                  cividis = viridis::cividis(n))

  ## 4) 连续型自定义渐变
  } else if (cmap %in% custom_maps) {
    if (cmap == "coolwarm") pal <- colorRampPalette(c("#3B4CC0","#FFFFFF","#B40426"))(n)
    if (cmap == "bwr")      pal <- colorRampPalette(c("#0000FF","#FFFFFF","#FF0000"))(n)
    if (cmap == "seismic")  pal <- colorRampPalette(c("#67001F","#FFFFFF","#053061"))(n)

  ## 5) 连续型 scico
  } else if (cmap %in% scico_maps) {
    pal <- scico::scico(n, palette = cmap)

  } else {
    stop("Unsupported colormap: ", cmap)
  }

  if (reverse) pal <- rev(pal)
  return(pal)
}
