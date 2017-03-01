
pheno.additive <- function(x, beta, sparse = TRUE) {
  a <- .Call("ps_vector_product",  PACKAGE='phenosim', x@bed, x@p, beta, sparse, FALSE)
  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0)
      names(a) <- x@ped$id
  }
  a
}

pheno.dominant <- function(x, beta, sparse = TRUE) {
  a <- .Call("ps_vector_product",  PACKAGE='phenosim', x@bed, x@p, beta, sparse, TRUE)
  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0)
      names(a) <- x@ped$id
  }
  a
}
  
