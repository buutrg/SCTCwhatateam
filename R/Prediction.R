#make sure to install dbscan: Density Based Clustering of Applications with Noise (DBSCAN) and Related Algorithms
# library(dbscan)

#input data:
#binary data for locations 3039*84
#binary data for cells 1297*84
#geometry data 3039*3
#gene list 60/40/20
#gdist

#output data:
#Top 10 posible locations of cells 1297*10

mcc <- function (actual, predicted)
{
  # Compute the Matthews correlation coefficient (MCC) score
  # Added zero denominator handling.
  # Avoided overflow error on large-ish products in denominator.
  #
  # actual = vector of true outcomes, 1 = Positive, 0 = Negative
  # predicted = vector of predicted outcomes, 1 = Positive, 0 = Negative
  # function returns MCC

  TP <- sum(actual == 1 & predicted == 1)
  TN <- sum(actual == 0 & predicted == 0)
  FP <- sum(actual == 0 & predicted == 1)
  FN <- sum(actual == 1 & predicted == 0)
  #TP;TN;FP;FN # for debugging
  sum1 <- TP+FP; sum2 <-TP+FN ; sum3 <-TN+FP ; sum4 <- TN+FN;
  denom <- as.double(sum1)*sum2*sum3*sum4 # as.double to avoid overflow error on large products
  if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {
    denom <- 1
  }
  mcc <- ((TP*TN)-(FP*FN)) / sqrt(denom)
  return(mcc)
}

# ranking locations and number of top 1 mcc
ranknmccMatirx <- function(bicells = NULL, biloc = NULL) {

  result <- NULL
  nmcc <- c()
  for (i in 1:nrow(bicells)) {
    mrow <- c()
    for (j in 1:nrow(biloc)) {
      preds <- bicells[i,]
      actuals <- biloc[j,]
      mrow <- c(mrow, mcc(preds, actuals))
    }
    rankl = order(mrow,decreasing=TRUE)
    nmcc <- c(nmcc, length(which(mrow == max(mrow))))
    result <- rbind(result, rankl)
  }
  #print ("OK")
  return(list(nmcc=nmcc,allrankl=result))
}

# search the neast neighbor to the selected list
nearbyl <- function(clist = NULL, exlist = NULL,geometry = NULL) {
  gdist <- as.matrix(dist(geometry))
  temp <- gdist[clist,exlist]
  index <- which(temp == min (temp), arr.ind = TRUE)[2]
  return (exlist[index])
}


# choose top 2 locations, then chose the closed 8 locations to them
top2mcc <- function(rankl = NULL,geometry = NULL) {
  result <- NULL
  for (i in 1:nrow(rankl)) {
    #  for (i in 3:3) {
    mrow <- rankl[i,1:2]
    for (j in 3:10) {
      nearl <- nearbyl(mrow,setdiff(rankl[i,],mrow),geometry)
      mrow <- c(mrow, nearl)
      #      print(mrow)
    }

    result <- rbind(result, mrow)
  }
  print ("OK")
  return(result)
}

densityl <- function(rankli = NULL, nmi = NULL,geometry = NULL){
  temp <- geometry[rankli[1:nmi],]
  scores <- lof(temp, k = 10)
  scorel <- order(scores,decreasing=FALSE)
  #print (rankli[scorel[1]])
  return (rankli[scorel[1:2]])
}


# choose locations based on density score
densmcc <- function(rankl = NULL,nmmcc = NULL, geometry = NULL) {

  result <- NULL
  for (i in 1:nrow(rankl)) {
    #   for (i in 3:3) {
    if (nmmcc[i] <= 10) {
      mrow <- rankl[i,1:2]
      for (j in 3:10) {
        nearl <- nearbyl(mrow,setdiff(rankl[i,],mrow),geometry)
        mrow <- c(mrow, nearl)
        #      print(mrow)
      }
    }
    else {
      mrow <- densityl(rankl[i,], nmmcc[i],geometry)
      for (j in 3:10) {
        nearl <- nearbyl(mrow,setdiff(rankl[i,],mrow),geometry)
        mrow <- c(mrow, nearl)
        #print(mrow)
      }
    }
    result <- rbind(result, mrow)
  }
  #print ("OK")
  return(result)
}
# calculate the Euclidist distance between vec1 and vec2
Euclidist <- function(vec1 = NULL, vec2 = NULL) {
  result <- sqrt((vec1[1]-vec2[1])^2+(vec1[2]-vec2[2])^2+(vec1[3]-vec2[3])^2)
  return (result)
}

#' QEvaluation
#'
#' @param realLoc real location.
#' @param predloc predicted locations.
#' @param geometry geometry data of locations.
#' @param plotDist TRUE if user wants to plot distribution of detailed results
#' @return sum of average distances.
#' @export
QEvaluation <- function(realLoc = NULL, predloc = NULL, geometry = NULL, plotDist = FALSE) {

  gdist <- as.matrix(dist(geometry))
  sum <- 0
  adists = NULL
  for (i in 1:nrow(realLoc)) {
    avgLoc <- 0
    for (j in 1:10){
      # calculate the distance between two locations
      d <- gdist[realLoc[i,1],predloc[i,j]]
      avgLoc <- avgLoc + d
    }
    avgLoc = avgLoc / 10
    adists = c(adists, avgLoc)
    sum = sum + avgLoc
  }
  if (plotDist)
    hist(adists)
  return(sum)
}
# library("ggpubr")

#
#' Listbycor
#'
#' Choose top 10 locations according cor
#'
#' @param normcells gene expression profiles in single cells.
#' @param normloc gene expression profiles in real locations.
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
#' @return 10 most likely locations.
#' @export
## sample * genes
# choose top 10 locations according cor
Listbycor <- function(normcells = NULL, normloc = NULL, method = NULL) {

  result <- NULL
  for (i in 1:nrow(normcells)) {
    mrow <- c()
    for (j in 1:nrow(normloc)) {
      dloc <- normcells[i,]
      dcell <- normloc[j,]
      mrow <- c(mrow,  cor(dloc, dcell,  method = method))
    }
    rankl = order(mrow,decreasing=TRUE)[1:10]
    result <- rbind(result, rankl)
  }
  return(result)
}
#
#' ListbyMCC
#'
#' Choose top 10 locations according MCC
#'
#' @param bicells84 binary gene expression profiles in single cells.
#' @param biloc84 binary gene expression profiles in real locations.
#' @param geneList driver genes.
#' @return 10 most likely locations.
#' @export
#input: binary data for cells, binary data for locations, gene list
#output: top 10 posible locations for each cell 1297*10
ListbyMCC <- function(bicells84 = NULL, biloc84 = NULL, geneList=NULL) {
  RankList <- ranknmccMatirx(bicells84[,geneList], biloc84[,geneList])
  result <- RankList$allrankl[,1:10]
  return (result)
}
#
#' ListbyMCCLOF
#'
#' Choose top 10 locations according MCC and LOF
#'
#' @param bicells84 binary gene expression profiles in single cells.
#' @param biloc84 binary gene expression profiles in real locations.
#' @param geometry coordinate location
#' @param geneList driver genes.
#' @return 10 most likely locations.
#' @export
#input: binary data for cells, binary data for locations, gene list
#output: top 10 posible locations for each cell 1297*10
ListbyMCCLOF <- function(bicells84 = NULL, biloc84 = NULL, geometry = NULL, geneList=NULL) {
  RankList <- ranknmccMatirx(bicells84[,geneList], biloc84[,geneList])
  result <- densmcc(RankList$allrankl, RankList$nmcc, geometry)
  return (result)
}
#' numofbins
#'
#' @param biExprs84 binary gene expression profiles in single cells or locations.
#' @param geneList driver genes.
#' @return numbers of unique bins in locations and cells.
#' @export
numofbins <- function(biExprs84 = NULL, geneList = NULL) {
  res <- length(table(do.call(paste, data.frame(biExprs84[,geneList]))))
  res
}
