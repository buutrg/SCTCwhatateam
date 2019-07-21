#' @import stats
#' @import graphics
#' @import reticulate
#' @import dbscan
#' @import geometry
#' @import ggplot2
#' @import dbscan
#' @import ggpubr
#' @import rgl
#' @import grDevices
#' @import CancerSubtypes
#' @import Rmagic

#'
# library(CancerSubtypes)

# Set environment variables here if any
# rootDir <- "C:/Users/phavy022/MyDoc/00Others/10Challenge"
# NoOfGenes <- 84
# Dorsal_ectoderm <- c("Ance", "CG2162", "Doc1", "Doc2", "egr", "peb", "tok", "ush", "zen")
# Neurectoderm <- c("ac", "brk", "CG8312", "l(1)sc", "mfas", "Ptp4E", "sog", "SoxN", "vnd")
# Mesoderm <- c("CG9005", "Cyp310a1", "GEFmeso", "ltl", "Mdr49", "Mes2", "NetA", "ry", "sna", "stumps", "twi", "wgn", "zfh1")
# Yolk_cells <- c("beat-IIIc", "CG8129", "CG8195", "Corp", "CNT1", "sisA", "ZnT77C")
# Pole_cells <- c("Pgc")

# Code from other files if any

#================================================================
# Functions if any
#================================================================

#'
processData <- function(dremi, Eucl) {

  rownames(dremi) <- rownames(Eucl)
  colnames(dremi) <- colnames(Eucl)
  dremi <- data.matrix(dremi)

  for (i in 1:(nrow(dremi)-1)) {
    for (j in (i+1):ncol(dremi)) {
      dremi[i,j] <- max(dremi[i,j], dremi[j,i])
      dremi[j, i] <- dremi[i,j]
    }
  }

  return(dremi)
}


#'
findGeneFarFromList <- function(l, Eucl) {
  r = NULL

  temp <- Eucl[l, !(colnames(Eucl) %in% l)]
  noOfRow <- length(row.names(temp))
  temp <- rbind(temp, 0)
  for (i in 1:(length(colnames(temp)))) {
    temp[noOfRow + 1, i] <- min(temp[1:noOfRow, i])
  }
  r <- colnames(temp)[which.max(temp[noOfRow + 1,])]

  return(r)
}

#'
findGeneDifFromList <- function(l, dremi) {
  r = NULL

  temp <- dremi[l, !(colnames(dremi) %in% l)]
  noOfRow <- length(row.names(temp))
  temp <- rbind(temp, 0)
  for (i in 1:(length(colnames(temp)))) {
    temp[noOfRow + 1, i] <- max(temp[1:noOfRow, i])
  }
  r <- colnames(temp)[which.min(temp[noOfRow + 1,])]

  return(r)
}



# @export
#'
getGeneList_M1 <- function(Eucl, sortedGeneDistanceList, expectedNo) {
  gl = NULL

  # Select 2 genes which are farest
  gl <- rbind(gl, sortedGeneDistanceList[1, 1])
  gl <- rbind(gl, sortedGeneDistanceList[1, 2])

  # Select the left
  while (length(gl) < expectedNo) {
    gl <- rbind(gl, findGeneFarFromList(gl, Eucl))
  }

  return(gl)
}

#'
getGeneList_M2 <- function(Eucl, seedList, expectedNo) {
  gl = NULL
  for (i in 1:length(seedList)) {
    gl <- rbind(gl, seedList[i])
  }

  # Select the left
  while (length(gl) < expectedNo) {
    gl <- rbind(gl, findGeneFarFromList(gl, Eucl))
  }

  return(gl)
}

#'
getGeneList_M61 <- function(dremi, sortedGeneRelationList, expectedNo) {
  gl = NULL

  # Select 2 genes which are not related most
  gl <- rbind(gl, sortedGeneRelationList[1, 1])
  gl <- rbind(gl, sortedGeneRelationList[1, 2])

  # Select the left
  while (length(gl) < expectedNo) {
    gl <- rbind(gl, findGeneDifFromList(gl, dremi))
  }

  return(gl)
}

#'
getGeneList_M62 <- function(dremi, seedList, expectedNo) {
  gl = NULL
  for (i in 1:length(seedList)) {
    gl <- rbind(gl, seedList[i])
  }

  # Select the left
  while (length(gl) < expectedNo) {
    gl <- rbind(gl, findGeneDifFromList(gl, dremi))
  }

  return(gl)
}
#' selectFarGenes
#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
#================================================================
# 1. Filter out similar genes
#================================================================

# setwd(rootDir)

# load(paste(rootDir, "/Data_names_fixed.RData", sep = ""))
# load(paste(rootDir, "/Euclidean.RData", sep = ""))
# View(head(Eucl))

# Build the list of gene distances

selectFarGenes = function(Eucl) {

  NoOfGenes = 84
  geneDistanceList = NULL
  for (r in 1:(NoOfGenes-1)) {
    for (c in (r+1):NoOfGenes) {
      geneDistanceList <- rbind(geneDistanceList, c(row.names(Eucl)[r], colnames(Eucl)[c], Eucl[r, c]))
    }
  }
  # View(head(geneDistanceList))
  # Sort
  sortedGeneDistanceList <- geneDistanceList[order(geneDistanceList[, 3], decreasing = TRUE), ]
  # View(head(sortedGeneDistanceList))

  # Select genes
  geneList20 <- getGeneList_M1(Eucl, sortedGeneDistanceList, 20)
  geneList40 <- getGeneList_M1(Eucl, sortedGeneDistanceList, 40)
  geneList60 <- getGeneList_M1(Eucl, sortedGeneDistanceList, 60)

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}

#' selectFarGenesWithSeeds
#' @param seedList Seed list
#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
#================================================================
# 2. Filter out similar genes with the seed in the groups
#================================================================

# setwd(rootDir)
#
# load(paste(rootDir, "/Euclidean.RData", sep = ""))
# View(head(Eucl))

# Get the seedList
selectFarGenesWithSeeds = function(seedList, Eucl){
  # seedList <- c(Dorsal_ectoderm, Neurectoderm, Mesoderm, Yolk_cells, Pole_cells)
  m <- match(seedList, colnames(Eucl))
  i <- which(m[] != "NA")
  seedList <- colnames(Eucl)[m[i]]

  # Select genes
  geneList20 <- getGeneList_M2(Eucl, seedList, 20)
  geneList40 <- getGeneList_M2(Eucl, seedList, 40)
  geneList60 <- getGeneList_M2(Eucl, seedList, 60)

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}




#================================================================
# 3. Use gene ranking
#================================================================
#' selectHighRakingGenes
#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
selectHighRakingGenes = function(Eucl){
  # setwd(rootDir)
  #
  # load(paste(rootDir, "/Euclidean.RData", sep = ""))
  # View(head(Eucl))

  a <- mean(Eucl)
  network <- Eucl
  for (i in 1:nrow(Eucl)) {
    for (j in 1:ncol(Eucl)) {
      if (i != j) {
        if (Eucl[i,j] <= a) {
          network[i,j] <- 1
        } else {
          network[i,j] <- 0
        }
      }
    }
  }

  # This code is adapted from the paper "Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and Expression Data"
  # Xu T, Le TD, Liu L, Wang R, Sun B, et al.
  # (2016) Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and Expression Data.
  # PLOS ONE 11(4): e0152792.
  # https://doi.org/10.1371/journal.pone.0152792
  colsum=colSums(network)
  colsum[which(colsum==0)]=1
  dev=t(apply(network, 1, function(x) x/colsum))
  num=nrow(dev)
  damping=0.85
  A=diag(num)-damping*dev
  B=solve(A)
  e=rep(1,num)
  delta=(1-damping)/num
  ranking=B %*% e * delta
  ranking=ranking/sum(ranking)

  # Order
  orderedRanking <- ranking[order(ranking[,1], decreasing = TRUE),]

  # Select genes
  geneList20 <- names(orderedRanking)[1:20]
  geneList40 <- names(orderedRanking)[1:40]
  geneList60 <- names(orderedRanking)[1:60]

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}




#================================================================
#' selectGenesByMADExpression
#  FSbyMAD with genes, cells
#================================================================
#' @param dge_normalized normalized single-cell matrix, rows for genes and columns for cells
#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
selectGenesByMADExpression = function(dge_normalized, Eucl){
  # setwd(rootDir)
  #
  # load(paste(rootDir, "/Euclidean.RData", sep = ""))
  # View(head(Eucl))
  # load(paste(rootDir, "/Data_names_fixed.RData", sep = ""))
  # View(head(dge_normalized))
  # Get 84 genes
  d <- dge_normalized[rownames(dge_normalized) %in% rownames(Eucl),]

  # Select genes
  geneList20 <- rownames(FSbyMAD(d, "topk", 20))
  geneList40 <- rownames(FSbyMAD(d, "topk", 40))
  geneList60 <- rownames(FSbyMAD(d, "topk", 60))

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))

}





#================================================================
#' selectGenesByMADDistance
#' Feature selection by MAD
# FSbyMAD with gene distances
#================================================================

#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
selectGenesByMADDistance = function(Eucl) {
  # setwd(rootDir)
  #
  # load(paste(rootDir, "/Euclidean.RData", sep = ""))
  # View(head(Eucl))

  # Select genes
  geneList20 <- rownames(FSbyMAD(Eucl, "topk", 20))
  geneList40 <- rownames(FSbyMAD(Eucl, "topk", 40))
  geneList60 <- rownames(FSbyMAD(Eucl, "topk", 60))

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}




#================================================================
#' selectGenesByInfluence
# 6. dremi
#================================================================
#' @param dremi distance matrix between 84 genes preprocessed by DREMI
#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'

selectGenesByInfluence = function(dremi, Eucl) {
  #================================================================
  # 6.1. dremi - Filter out similar genes
  #================================================================

  # setwd(rootDir)
  #
  # load(paste(rootDir, "/Euclidean.RData", sep = ""))
  # View(head(Eucl))

  # dremi <- read.csv(file="./dremi1.csv", header=TRUE, row.names=1)

  NoOfGenes = nrow(Eucl)

  dremi <- processData(dremi, Eucl)

  # Build the list of gene relations
  geneRelationList = NULL
  for (r in 1:(NoOfGenes-1)) {
    for (c in (r+1):NoOfGenes) {
      geneRelationList <- rbind(geneRelationList, c(row.names(dremi)[r], colnames(dremi)[c], dremi[r, c]))
    }
  }
  # View(head(geneRelationList))
  # Sort
  sortedGeneRelationList <- geneRelationList[order(geneRelationList[, 3], decreasing = FALSE), ]
  # View(head(sortedGeneRelationList))

  # Select genes
  geneList20 <- getGeneList_M61(dremi, sortedGeneRelationList, 20)
  geneList40 <- getGeneList_M61(dremi, sortedGeneRelationList, 40)
  geneList60 <- getGeneList_M61(dremi, sortedGeneRelationList, 60)

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}



#================================================================
#' selectGenesByInfluenceWithSeeds
# 6.2. dremi - Filter out similar genes with the seed in the groups
#================================================================
#' @param seedList the initial seed list of genes
#' @param dremi distance matrix between 84 genes preprocessed by DREMI
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'

selectGenesByInfluenceWithSeeds = function(seedList, dremi){
  # setwd(rootDir)
  #
  # load(paste(rootDir, "/Euclidean.RData", sep = ""))
  # View(head(Eucl))

  # dremi <- read.csv(file="./dremi.csv", header=TRUE, row.names=1)
  # dremi <- processData(dremi, Eucl)

  # Get the seedList
  # seedList <- c(Dorsal_ectoderm, Neurectoderm, Mesoderm, Yolk_cells, Pole_cells)
  m <- match(seedList, colnames(dremi))
  i <- which(m[] != "NA")
  seedList <- colnames(dremi)[m[i]]

  # Select genes
  geneList20 <- getGeneList_M62(dremi, seedList, 20)
  geneList40 <- getGeneList_M62(dremi, seedList, 40)
  geneList60 <- getGeneList_M62(dremi, seedList, 60)

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}



#================================================================
#' selectHighRakingGenesByInfluence
# 6.3. dremi - Use gene ranking
#================================================================
#' @param dremi distance matrix between 84 genes preprocessed by DREMI
#' @param Eucl Euclidean distance matrix between 84 genes
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
selectHighRakingGenesByInfluence = function(dremi, Eucl) {
  # load(paste(rootDir, "/Euclidean.RData", sep = ""))
  # View(head(Eucl))

  # dremi <- read.csv(file="./dremi.csv", header=TRUE, row.names=1)
  dremi <- processData(dremi, Eucl)

  a <- mean(dremi)
  network <- dremi
  for (i in 1:nrow(dremi)) {
    for (j in 1:ncol(dremi)) {
      if (i != j) {
        if (dremi[i,j] >= a) {
          network[i,j] <- 1
        } else {
          network[i,j] <- 0
        }
      }
    }
  }

  # This code is adapted from the paper "Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and Expression Data"
  # Xu T, Le TD, Liu L, Wang R, Sun B, et al.
  # (2016) Identifying Cancer Subtypes from miRNA-TF-mRNA Regulatory Networks and Expression Data.
  # PLOS ONE 11(4): e0152792.
  # https://doi.org/10.1371/journal.pone.0152792
  colsum=colSums(network)
  colsum[which(colsum==0)]=1
  dev=t(apply(network, 1, function(x) x/colsum))
  num=nrow(dev)
  damping=0.85
  A=diag(num)-damping*dev
  B=solve(A)
  e=rep(1,num)
  delta=(1-damping)/num
  ranking=B %*% e * delta
  ranking=ranking/sum(ranking)

  # Order
  orderedRanking <- ranking[order(ranking[,1], decreasing = TRUE),]

  # Select genes
  geneList20 <- names(orderedRanking)[1:20]
  geneList40 <- names(orderedRanking)[1:40]
  geneList60 <- names(orderedRanking)[1:60]

  # Write files
  return(list("ls20"=geneList20, "ls40"=geneList40, "ls60"=geneList60))
}


