
#' preproc_MAGIC
#'
#' Preprocessing expression data by MAGIC and get distance by Euclidean
#'
#' @param data expression data with row as genes columns as samples or cells
#' @param gene_markers A list of gene markers e.g. 84 genes, to get distances
#' @return Distance between gene markers
#' @export
#'
preproc_MAGIC = function(data, gene_markers) {
  after_magic = magic(data, gene_markers, t = 'auto', seed = 1234)
  # mat_magic = as.matrix(after_magic)
  # data_new = mat_magic[gene_markers, ]
  Eucl = as.matrix(dist(t(after_magic$result)))
  return(Eucl)
}

#' preproc_MAGIC_dremi
#'
#' Preprocessing expression data by MAGIC and get distance by dremi
#'
#' @param data Expression data with row as genes columns as samples or cells
#' @param gene_markers A list of gene markers e.g. 84 genes, to get distances
#' @return Distance between gene markers
#' @export
#'
preproc_MAGIC_dremi = function(data, gene_markers) {

  scprep = import("scprep")
  after_magic = magic(t(data), t = 'auto', k=10, npca=420, seed = 1234)
  mat_magic = as.matrix(after_magic)
  data_new = mat_magic[gene_markers, ]

  distance = apply(data_new, 1, function(x) {
    apply(data_new, 1, function(y) {
      scprep$stats$knnDREMI(x, y)
    })
  })

  return(distance)
}

