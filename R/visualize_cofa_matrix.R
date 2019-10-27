##----------------------------------------------------------------------------------------------
## Functions for Visualization
##----------------------------------------------------------------------------------------------

# parameters:
#   mat - matrix object
# returns:
#   lowere triangular version of matrix
#
get_lower_tri <- function(mat){
  mat[upper.tri(mat)] <- NA
  return(mat)
}

# parameters:
#   mat - matrix object
# returns:
#   matrix with labels ordered so that similar values are grouped together
reorder_mat <- function(freqMat){
  hc <- cluster_mat(freqMat)
  return(freqMat[hc$order, hc$order])
}

# parameters:
#   mat - matrix object
#   verbose - if true, print distance matrix
# returns:
#   returns hclust object
cluster_mat <- function(freqMat, verbose=FALSE, na.value=0){
  freqMat[is.na(freqMat)] <- na.value
  dd <- as.dist((1-freqMat)/2)
  if ((verbose) == TRUE){print(dd)} # for debugging
  return(stats::hclust(dd))
}

# summary:melt a matrix and order the factors so they will be plotted in the
# same order as they are given
# parameters:
#   temp - matrix that is ordered or unordered
# returns: melted matrix ready for ggplot
meltForViz <- function(temp){

  # melt it
  melted_freq_mat <- melt(temp,
                          na.rm = TRUE,
                          varnames = c("var1","var2"),
                          value.name ="value")

  # make the factors ordered so they will plot in right order
  melted_freq_mat$var1 <- factor(melted_freq_mat$var1, levels = c(colnames(temp)) )
  melted_freq_mat$var2 <- factor(melted_freq_mat$var2, levels = c(colnames(temp)) )

  return(melted_freq_mat)
}


# summary: Creates ggplot visualization of co-frequency matrix
# parameters:
#   a frequency matrix object with atleast the columns labelled
# returns:
#   ggplot object of cofrequency matrix plot
#
vizCoFreqMat <- function(freqMat, order=TRUE, alph=FALSE, sigfigs=1, text=TRUE){

  if ( sum(!is.na(freqMat)) == 0 ) stop("given frequency matrix is all NA")

  # sorting rows/columns alphabetically
  if (alph) {
    freqMat <- freqMat[rev(sort(rownames(freqMat))),rev(sort(rownames(freqMat)))]
  }

  # reformat frequency matrix
  if (order){
    temp <- get_lower_tri(reorder_mat(freqMat))
  } else {
    temp <- get_lower_tri(freqMat)
  }

  melted_freq_mat = meltForViz(temp)

  # plot with ggplot
  p <- ggplot(data = melted_freq_mat, aes(var1, var2, fill = value, label=signif(value,sigfigs))) +
    geom_tile(color="white") +
    coord_fixed() + theme_light() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = c(0.4, 0.9),
          legend.direction = "horizontal",
          axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0)) +
    guides(fill = guide_colorbar(barwidth = 10, barheight = 0.5)) +
    scale_fill_viridis(name="", direction=1,
                       begin=0, end=0.95, limits=c(0,1), alpha=0.5)

  if (text) { p <- p + geom_text(size=1, alpha=0.7) }

  return(p)
}
