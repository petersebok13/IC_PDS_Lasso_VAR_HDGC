######### Packages dependencies:  ##########
library(igraph)
library(zoo)
library(qpcR)
library(dplyr)
library(glmnet)
library(aod)
library(readxl)
############################################

NetGC_HDGC <- function(data, log = FALSE, std = FALSE,
                       alpha = 1, sign = 0.05, crit = "bic",
                       lags = 1, method = c("ic", "theoretical", "cv"),
                       plot = c("circle", "veins"), cluster = FALSE, verbose = TRUE,
                       title = "HD Granger Volatility Network") {
  
  method <- match.arg(method)
  if (!is.data.frame(data)) stop("data must be a dataframe")
  if (log) data <- log(data)
  if (std) data <- as.data.frame(scale(data))
  
  N <- ncol(data)
  varnames <- colnames(data)
  adj_matrix <- matrix(0, nrow = N, ncol = N)
  rownames(adj_matrix) <- varnames
  colnames(adj_matrix) <- varnames
  
  for (i in 1:N) {
    for (j in 1:N) {
      if (i == j) next
      
      # Apply HDGC function
      pval <- tryCatch({
        HDGC_x_y(x = varnames[j], y = varnames[i], data = data,
             std = FALSE, crit = crit, alpha = alpha,
             sign = sign, method = method)
      }, error = function(e) NA)
      
      if (!is.na(pval) && pval < sign) {
        adj_matrix[j, i] <- 1
        if (verbose) cat(paste0(varnames[j], " → ", varnames[i], " (p = ", signif(pval, 4), ")\n"))
      } else {
        if (verbose) cat(paste0(varnames[j], " -/→ ", varnames[i], " (p = ", signif(pval, 4), ")\n"))
      }
    }
  }
  
  # Build graph
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", diag = FALSE)
  V(g)$label <- V(g)$name
  
  if (plot[1] == "circle") {
    plot(g, layout = layout_in_circle, main = title,
         edge.arrow.size = .1, vertex.color = "gold", vertex.size = 5,
         vertex.frame.color = "gray", vertex.label.color = "black",
         vertex.label.cex = 1.5, vertex.label.dist = 0.75, edge.curved = 0)
  } else if (plot[1] == "veins") {
    plot(g, layout = layout_nicely, vertex.size = 17.5,
         edge.arrow.size = .15, vertex.label.cex = 0.70,
         vertex.label.dist = 0, main = title)
  }
  
  if (cluster) {
    g_undirected <- as.undirected(g, mode = "collapse")
    ceb <- cluster_edge_betweenness(g_undirected, directed = TRUE)
    dendPlot(ceb, mode = "hclust")
    plot(ceb, g_undirected, vertex.size = 25,
         vertex.label.color = "black", vertex.label.cex = 0.5, edge.curved = 0)
  }
  
  invisible(adj_matrix)
}
