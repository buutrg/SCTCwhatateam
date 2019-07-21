# library(rgl)

myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

# 3D plot positions of predicted locations.
#' plot3DCellPositions
#'
#' @param geometry geometry data of locations.
#' @param locidx Numeric vector of predicted locations.
#' @export
plot3DCellPositions <- function(geometry = NULL, locidx = NULL) {
  #locidx <- sample(1:3039, 1297, replace=T)
  freqidx = as.data.frame(table(locidx))
  locidx = as.numeric(levels(freqidx$locidx))[freqidx$locidx]
  predloc = geometry[locidx,]
  cols <- myColorRamp(c("grey", "black"), freqidx$Freq)
  plot3d(predloc,xlab = "xcoord", ylab = "ycoord", zlab = "zcoord", col = cols)
}

# 2D plot positions of predicted locations.
#' plot2DCellPositions
#'
#' @param geometry geometry data of locations.
#' @param locidx Numeric vector of predicted locations.
#' @param xcoord x coordinate
#' @param zcoord z coordinate
#' @export
plot2DCellPositions <- function(geometry = NULL, locidx = NULL, xcoord, zcoord) {
  #locidx <- sample(1:3039, 1297, replace=T)
  freqidx = as.data.frame(table(locidx))
  locidx = as.numeric(levels(freqidx$locidx))[freqidx$locidx]
  predloc = geometry[locidx,]
  # if (require(ggplot2)) {
  ggplot(predloc) +
    geom_point(aes(x=xcoord, y=zcoord, color=freqidx$Freq)) +
    labs(color='Frequence')
  # }
}

# 3D plot gene patterns of predicted locations.
#' plot3DgenePattern
#'
#' @param geometry geometry data of locations.
#' @param locidx Numeric vector of predicted locations.
#' @param geneProfile gene profile of cells
#' @export
plot3DgenePattern <- function(geometry = NULL, locidx = NULL, geneProfile = NULL) {
  #locidx <- sample(1:3039, 1297, replace=T)
  temp = t(rbind(locidx, geneProfile))
  rownames(temp) = temp[,1]
  temp =apply(temp, 2, function(x) tapply(x, rownames(temp), mean))
  predloc = geometry[temp[,1],]
  cols <- myColorRamp(c("yellow", "red"), temp[,2])
  plot3d(predloc,xlab = "xcoord", ylab = "ycoord", zlab = "zcoord", col = cols)
}

# 2D plot gene patterns of predicted locations.
#' plot2DgenePattern
#'
#' @param geometry geometry data of locations.
#' @param locidx Numeric vector of predicted locations.
#' @param geneProfile gene profile of cells
#' @param xcoord x coordinate
#' @param zcoord z coordinate
#' @export
plot2DgenePattern <- function(geometry = NULL, locidx = NULL, geneProfile = NULL, xcoord, zcoord) {
  #locidx <- sample(1:3039, 1297, replace=T)
  temp = t(rbind(locidx, geneProfile))
  rownames(temp) = temp[,1]
  temp =apply(temp, 2, function(x) tapply(x, rownames(temp), mean))
  predloc = geometry[temp[,1],]
  # if (require(ggplot2)) {
  ggplot(predloc) +
    geom_point(aes(x=xcoord, y=zcoord, color=temp[,2])) +
    labs(color='gene expression')
  # }
}
