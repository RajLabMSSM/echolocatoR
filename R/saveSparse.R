#' Save LD matrix as a sparse matrix
#'
#' Converting LD matrices to sparse format reduces file size by half.
#' @family LD
#' @keywords internal
saveSparse <- function(LD_matrix,
                       LD_path,
                       verbose=T){
  # https://cmdlinetips.com/2019/05/introduction-to-sparse-matrices-in-r/
  LD_sparse <- Matrix::Matrix(as.matrix(LD_matrix), sparse = T)
  # printer("Dense size:")
  # print(object.size(LD_matrix),units="auto")
  # printer("Sparse size:")
  # print(object.size(LD_sparse),units="auto")
  printer("LD:: Saving LD as sparse matrix ==>",LD_path,v=verbose)
  saveRDS(LD_sparse, LD_path)
}
