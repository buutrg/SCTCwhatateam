# bdtnp_file = "/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/DC_package/data/bdtnp.csv"
# dge_normalized_file = "/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/DC_package/data/dge_normalized.csv"
#
#
# linGen = function(bdtnp_file, dge_normalized_file)

#' build_lin
#' @param lin_dir directory to python code
#' @return 3 lists of 20, 40, 60 selected genes
#' @export
#'
#'
build_lin = function(lin_dir) {
  # library("reticulate")
  # lin_dir = "/Users/buutruong/TMTB/Computational_Research/UniSA/Dream_Challenges/SCTC/DC_package/package/SCTCwhatateam/linMethods_py/linmethods.py"
  source_python(lin_dir)
}
