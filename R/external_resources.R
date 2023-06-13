#########################   External Resources Module   ########################
##
## @description
## This module contains functions adapted from methods described in the
## following paper:
##          (https://doi.org/10.1016/j.celrep.2016.10.057).
##   - SCDC: Dong, Meichen, et al. "SCDC: Bulk Gene Expression Deconvolution by
##           Multiple Single-Cell RNA Sequencing References." Briefings in
##           Bioinformatics, vol. 22, no. 1, 2021, pp. 416-427.
##           [DOI: 10.1093/bib/bbz166](https://doi.org/10.1093/bib/bbz166).
##
## These functions are a part of the Decoflex package. The original
## implementations of these functions can be found in the respective GitHub
## repositories as cited in the documentation of each function within this
## module.
##
## @note
## This module and its functions should be used in accordance with the
## respective licenses and guidelines provided by the authors of the original
## methods.
##

#' SCDC_basis_final
#'
#' @description
#' This function constructs a basis matrix based on an ExpressionSet object for
#' single cells. It is based on the function created in the SCDC R
#' implementation. Please refer to the original implementation for more
#' information: https://meichendong.github.io/SCDC/. Version: 0.0.0.9000
#'
#' @param x An ExpressionSet object for single cells. This object will serve as
#'  the basis for the construction of the basis matrix.
#' @param ct.sub A vector of cell types that are selected to construct the basis
#'  matrix. This parameter is optional, with the default being NULL.
#' @param ct.varname A string that indicates the variable name for 'cell types'
#'  in the 'x' dataset.
#' @param sample A string that indicates the variable name for subject/samples
#'  in the 'x' dataset.
#' @param ct.cell.size A vector of cell size factors corresponding to the ct.sub
#'  according to prior knowledge. The default is NULL, which means the "library
#'  size" is calculated based on the data. If a vector is provided, the names
#'  of the vector (names(ct.cell.size)) should not be NULL.
#' @param verbose A boolean flag for controlling the amount of information
#'  printed during the execution of the function. Default is FALSE.
#'
#' @return
#' A list with the following components:
#'   - Basis matrix
#'   - Sum of cell-type-specific library size
#'   - Sample variance matrix
#'   - Basis matrix by mvw
#'   - mvw matrix.
#'
#' @import Biobase
#' @import stats
#'
#' @references
#' SCDC: a proportion estimation method for single cell RNA-Seq data.
#' <https://meichendong.github.io/SCDC/>.
#'
#' @export
SCDC_basis_final <- function(x, ct.sub = NULL, ct.varname, sample,
                             ct.cell.size = NULL, verbose = FALSE){
  # Select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])
  }
  ct.sub <- ct.sub[!is.na(ct.sub)]
  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]

  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(Biobase::exprs(x.sub)) > 0,]

  # Calculate sample mean & sample variance matrix: genes by cell types
  countmat <- Biobase::exprs(x.sub)
  ct.id <- droplevels(as.factor(x.sub@phenoData@data[,ct.varname]))
  sample.id <- as.character(x.sub@phenoData@data[,sample])
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  sigma <- sapply(unique(mean.id[, 1]), function(id) {
    y = mean.mat[, mean.id[, 1] %in% id]
    if (is.null(dim(y))){
      res = rep(0, length(y))
      message("Warning: the cell type [", id,"] is only available",
              " in at most 1 subject!")
    } else {
      res = apply(y, 1, stats::var, na.rm = TRUE)
    }
    return(res)
  })

  sum.mat2 <- sapply(unique(sample.id), function(sid) {
    sapply(unique(ct.id), function(id) {
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in% sid])
      if (ncol(y)>0){
        out = sum(y)/ncol(y)
      } else {
        out = 0
      }
      return(out)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)

  # library size factor calculated from the samples:
  if (is.null(ct.cell.size)){
    sum.mat <- rowMeans(sum.mat2, na.rm = T)
  } else {
    if (is.null(names(ct.cell.size))){
      stop("Cell size factor vector requires cell type names...")
    } else {
      sum.mat <- ct.cell.size
    }
  }

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x, ...) {
    y <- apply(x, 1, max, na.rm = TRUE)
    if (stats::median(y, na.rm = T) == 0){
      outx = y
    }else{
      outx = y/stats::median(y, na.rm = T)
    }
    return(outx)
  }

  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      if (ncol(y)>0){
        out = apply(y, 1, stats::var, na.rm = T)
      } else {
        out = rep(0, nrow(y))
      }
      return(out)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  q15 <- apply(var.adj, 2, function(zz) {
    z1 = min(zz[zz > 0])
    z2 = stats::quantile(zz, 0.15, na.rm = T)
    return(max(z1, z2))
  })
  q85 <- apply(var.adj, 2, stats::quantile, probs = 0.85, na.rm = T)
  var.adj.q <- t(apply(var.adj, 1, function(y) {
    y[y < q15] <- q15[y < q15]
    y[y > q85] <- q85[y > q85]
    return(y)
  }))

  var.adj.q <- t(apply(var.adj, 1, function(y){
    y[y<q15] <- q15[y<q15]
    y[y>q85] <- q85[y>q85]
    return(y)}
  )
  )

  if(verbose){
    message("Creating Basis Matrix adjusted for maximal variance weight")
  }

  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    sid = unlist(strsplit(id,'%'))[2]
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q[,sid]), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,ct.sub]
  sigma <- sigma[, ct.sub]
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

#' SCDC_basis_ONE_final
#'
#' @description
#' This function constructs a basis matrix for single cells from one subject,
#' based on an ExpressionSet object. It is based on the function created in the
#' SCDC R implementation. Please refer to the original implementation for more
#' information: https://meichendong.github.io/SCDC/. Version: 0.0.0.9000
#'
#' @param x An ExpressionSet object for single cells. This object will serve as
#'  the basis for the construction of the basis matrix.
#' @param ct.sub A vector of cell types that are selected to construct the basis
#'  matrix. This parameter is optional, with the default being NULL.
#' @param ct.varname A string that indicates the variable name for 'cell types'
#'  in the 'x' dataset.
#' @param sample A string that indicates the variable name for subject/samples
#'  in the 'x' dataset.
#' @param ct.cell.size A vector of cell size factors corresponding to the ct.sub
#'  according to prior knowledge. The default is NULL, which means the "library
#'  size" is calculated based on the data. If a vector is provided, the names
#'  of the vector (names(ct.cell.size)) should not be NULL.
#' @param verbose A boolean flag for controlling the amount of information
#'  printed during the execution of the function. Default is FALSE.
#'
#' @return
#' A list with the following components:
#'   - Basis matrix
#'   - Sum of cell-type-specific library size
#'   - Sample variance matrix
#'   - Basis matrix by mvw
#'   - mvw matrix.
#'
#' @import Biobase
#' @import stats
#'
#' @references
#' SCDC: a proportion estimation method for single cell RNA-Seq data.
#' <https://meichendong.github.io/SCDC/>.
#'
#' @export
SCDC_basis_ONE_final <- function(x , ct.sub = NULL, ct.varname, sample,
                                 ct.cell.size = NULL, verbose = FALSE){
  # select only the subset of cell types of interest
  if (is.null(ct.sub)){
    ct.sub <- unique(x@phenoData@data[,ct.varname])[!is.na(unique(x@phenoData@data[,ct.varname]))]
  }else{}
  ct.sub <- ct.sub[!is.na(ct.sub)]

  x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]
  # qc: remove non-zero genes
  x.sub <- x.sub[rowSums(Biobase::exprs(x.sub)) > 0,]

  # Calculate sample mean & sample variance matrix: genes by cell types
  countmat <- Biobase::exprs(x.sub)

  ct.id <- x.sub@phenoData@data[,ct.varname]
  sample.id <- x.sub@phenoData@data[,sample]
  ct_sample.id <- paste(ct.id,sample.id, sep = '%')

  mean.mat <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    apply(y,1,sum, na.rm = TRUE)/sum(y)
  })
  mean.id <- do.call('rbind',strsplit(unique(ct_sample.id), split = '%'))

  # By subj, then take avg
  sum.mat2 <- sapply(unique(sample.id), function(sid) {
    sapply(unique(ct.id), function(id) {
      y = as.matrix(countmat[, ct.id %in% id & sample.id %in%
                               sid])
      if (ncol(y)>0){
        out = sum(y)/ncol(y)
      } else {
        out = 0
      }
      return(out)
    })
  })
  rownames(sum.mat2) <- unique(ct.id)
  colnames(sum.mat2) <- unique(sample.id)

  if (is.null(ct.cell.size)){
    sum.mat <- rowMeans(sum.mat2, na.rm = T)
  } else {
    if (is.null(names(ct.cell.size))){
      stop("Cell size factor vector requires cell type names.")
    } else {
      sum.mat <- ct.cell.size
    }
  }

  basis <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1, mean, na.rm = TRUE)
  })

  # weighted basis matrix
  my.max <- function(x, ...) {
    y <- apply(x, 1, max, na.rm = TRUE)
    if (stats::median(y, na.rm = T) == 0){
      outx = y
    }else{
      outx = y/stats::median(y, na.rm = T)
    }
    return(outx)
  }

  var.adj <- sapply(unique(sample.id), function(sid) {
    my.max(sapply(unique(ct.id), function(id) {
      y = countmat[, ct.id %in% id & sample.id %in% sid,
                   drop = FALSE]
      if (ncol(y)>0){
        out = apply(y, 1, stats::var, na.rm = T)
      } else {
        out = rep(0, nrow(y))
      }
      return(out)
    }), na.rm = T)
  })
  colnames(var.adj) <- unique(sample.id)

  q15 <- apply(var.adj, 2, function(zz){
    z1 = min(zz[zz>0])
    z2 = stats::quantile(zz, 0.15, na.rm = T)
    return(max(z1,z2))
  })
  q85 <- apply(var.adj,2, stats::quantile, probs = 0.85, na.rm =T)

  var.adj.q <- as.matrix(apply(var.adj, 1,
                               function(y){y[y<q15] <- q15[y<q15]
                               y[y>q85] <- q85[y>q85]
                               return(y)}) ) #+ 1e-4

  if(verbose){
    message("Creating Basis Matrix adjusted for maximal variance weight.")
  }

  mean.mat.mvw <- sapply(unique(ct_sample.id), function(id){
    y = as.matrix(countmat[, ct_sample.id %in% id])
    yy = sweep(y, 1, sqrt(var.adj.q), '/')
    apply(yy,1,sum, na.rm = TRUE)/sum(yy)
  })

  basis.mvw <- sapply(unique(mean.id[,1]), function(id){
    z <- sum.mat[mean.id[,1]]
    mean.mat.z <- t(t(mean.mat.mvw)*z)
    y = as.matrix(mean.mat.z[,mean.id[,1] %in% id])
    apply(y,1,mean, na.rm = TRUE)
  })

  # reorder columns
  basis.mvw <- basis.mvw[,levels(ct.sub)[ct.sub]]

  # In the one subject case, no variance is calculated.
  sigma <- NULL
  basis <- basis[, ct.sub]
  sum.mat <- sum.mat[ct.sub]

  return(list(basis = basis, sum.mat = sum.mat,
              sigma = sigma, basis.mvw = basis.mvw, var.adj.q = var.adj.q))
}

#' decoflex_build_cell_reference
#'
#' @description
#' This function constructs a reference matrix for single cells, based on an
#' ExpressionSet object. It is based on the function created in the SCDC R
#' implementation. Please refer to the original implementation for more
#' information: https://meichendong.github.io/SCDC/. Version: 0.0.0.9000
#'
#' @param x An ExpressionSet object for single cells. This object will serve as
#'  the basis for the construction of the reference matrix.
#' @param ct.sub A vector of cell types that are selected to construct the
#'  reference matrix. This parameter is optional, with the default being NULL.
#' @param ct.varname A string that indicates the variable name for 'cell types'
#'  in the 'x' dataset.
#' @param sample A string that indicates the variable name for subject/samples
#'  in the 'x' dataset.
#' @param ct.cell.size A vector of cell size factors corresponding to the
#'  ct.sub according to prior knowledge. The default is NULL, which means the
#'  "library size" is calculated based on the data. If a vector is provided,
#'  the names of the vector (names(ct.cell.size)) should not be NULL.
#' @param verbose A boolean flag for controlling the amount of information
#'  printed during the execution of the function. Default is FALSE.
#'
#' @return
#' A list with the following components:
#'   - basis: The constructed reference matrix.
#'   - detailed: The detailed output from SCDC_basis_ONE_final or
#'      SCDC_basis_final, which includes:
#'       - Basis matrix
#'       - Sum of cell-type-specific library size
#'       - Sample variance matrix
#'       - Basis matrix by mvw
#'       - mvw matrix.
#'
#' @references
#' SCDC: a proportion estimation method for single cell RNA-Seq data.
#' <https://meichendong.github.io/SCDC/>.
#'
#' @export
decoflex_build_cell_reference <- function(x, ct.sub = NULL, ct.varname,
                                          sample, ct.cell.size = NULL,
                                          verbose = FALSE){

  if(verbose){
    print('Parameters decoflex_build_cell_reference function:')
    print(paste0('ct.sub: ', ct.sub))
    print(paste0('ct.varname: ', ct.varname))
    print(paste0('sample: ', sample))
  }

  if(verbose){
    message('Number of cells for the level: ')

    # For each sample I calculate the number of cells presents
    for (counter in unique(x[[sample]])) {
      print(paste0('sample: ', counter))

      # In case that ct.sub is null, I will take all the values of the
      # variable sample
      if (is.null(ct.sub)){
        ct.sub <- unique(x@phenoData@data[,ct.varname])[!is.na(unique(x@phenoData@data[,ct.varname]))]
      }else{}

      # Calculation of the count function over the different samples
      summary.table <- table(x[[ct.varname]][x[[sample]]==counter & x[[ct.varname]] %in% ct.sub])
      summary.table <- summary.table[summary.table>0]
      print(summary.table)
    }
  }

  # Reference that is going to be returned.
  reference = NULL

  # Calculating the possible samples??
  use_single_sampling_reference = FALSE
  if(is.null(sample)){
    use_single_sampling_reference = TRUE
  }else{

    # Select only the subset of cell types of interest
    if (is.null(ct.sub)){
      ct.sub <- unique(x@phenoData@data[,ct.varname])[!is.na(unique(x@phenoData@data[,ct.varname]))]
    }else{}
    ct.sub <- ct.sub[!is.na(ct.sub)]

    x.sub <- x[,x@phenoData@data[,ct.varname] %in% ct.sub]

    # Check the possible values of the samples
    sample.id <- x.sub@phenoData@data[,sample]
    samples_values <- unique(sample.id)

    # Verify the possible samples
    if(length(samples_values) > 1){
      use_single_sampling_reference = FALSE
    }else{
      use_single_sampling_reference = TRUE
    }
  }

  # Check which analysis I should use.
  result_reference <- NULL
  if(use_single_sampling_reference){

    if(verbose){
      message(paste0("Creating reference without samples for: ",
                     paste(shQuote(ct.sub), collapse=", ")))
    }

    result_reference <- SCDC_basis_ONE_final(x = x, ct.sub = ct.sub,
                                             ct.varname = ct.varname,
                                             sample = sample,
                                             ct.cell.size = ct.cell.size,
                                             verbose = verbose)
    # In this case I have to use the basis matrix
    reference <- result_reference$basis
  }else{

    if(verbose){
      message(paste0("Creating reference with samples for: ",
                     paste(shQuote(ct.sub), collapse=", ")))
    }

    result_reference <- SCDC_basis_final(x = x, ct.sub = ct.sub,
                                         ct.varname = ct.varname,
                                         sample = sample,
                                         ct.cell.size = ct.cell.size,
                                         verbose = verbose)

    # In this case I have to use the basis.mvw that takes in account the
    # different samples
    reference <- result_reference$basis.mvw
  }

  # Return of the matrix with the reference of celltypes
  return(list(basis = reference, detailed = result_reference))
}

#' generateBulk_allcells
#'
#' @description
#' This function constructs Pseudo bulk samples for single cells from one
#' subject based on an ExpressionSet object. It is based on the function
#' created in the SCDC R implementation. Please refer to the original
#' implementation for more information: https://meichendong.github.io/SCDC/.
#' Version: 0.0.0.9000
#'
#' @param eset An ExpressionSet object for single cells. This object will serve
#'  as the basis for the construction of the Pseudo bulk samples.
#' @param ct.varname A string that indicates the variable name for 'cell types'
#'  in the 'eset' dataset.
#' @param sample A string that indicates the variable name for subject/samples
#'  in the 'eset' dataset.
#' @param disease A string indicating the health condition of subjects.
#' @param ct.sub A vector of cell types that are selected to construct Pseudo
#'  bulk samples. This parameter is optional, with the default being NULL. If
#'  NULL, then all cell types are used.
#' @param verbose A boolean flag for controlling the amount of information
#'  printed during the execution of the function. Default is FALSE.
#'
#' @return
#' A list with the following components:
#'   - Pseudo bulk samples ExpressionSet
#'   - Actual cell-type proportions
#'
#' @import Biobase
#' @import Matrix
#' @import stats
#'
#'
#' @references
#' SCDC: a proportion estimation method for single cell RNA-Seq data.
#' <https://meichendong.github.io/SCDC/>.
#'
#' @export
generateBulk_allcells <- function(eset, ct.varname, sample, disease = NULL,
                                  ct.sub = NULL, verbose = FALSE){

  if(verbose){
    print('Parameters generateBulk_allcells function:')
    print(paste0('ct.sub: ', ct.sub))
    print(paste0('ct.varname: ', ct.varname))
    print(paste0('sample: ', sample))
  }

  if (is.null(ct.sub)){
    ct.sub <- unique(eset@phenoData@data[,ct.varname])
  }
  eset <- eset[, eset@phenoData@data[,ct.varname] %in% ct.sub]
  cluster.id <- eset@phenoData@data[,ct.varname]
  sample.id <- eset@phenoData@data[,sample]
  condition.id <- eset@phenoData@data[,disease]

  # expression
  pseudo.exprs <- sapply(unique(sample.id), function(sid){
    y <- Biobase::exprs(eset)[, sample.id %in% sid]
    Matrix::rowSums(y, na.rm = T)
  })
  colnames(pseudo.exprs) <- unique(sample.id)
  # true proportion: sample by cell types
  ncount <- table(sample.id, cluster.id)
  true.prop <- ncount / Matrix::rowSums(ncount, na.rm = T)
  true.prop <- true.prop[stats::complete.cases(true.prop),]

  # eset for pseudo bulk sample
  if (is.null(disease)){
    pseudo.disease <- NA
  } else {
    pseudo.disease <- sapply(unique(sample.id), function(sid){
      condition.id[sample.id == sid][1]
    })
  }
  pseudo.pdata <- data.frame(sample = colnames(pseudo.exprs),
                             disease = pseudo.disease)
  pseudo.fdata <- data.frame(genes = rownames(pseudo.exprs))
  rownames(pseudo.fdata) <- rownames(pseudo.exprs)
  pseudo_eset <- getESET(exprs = pseudo.exprs,
                         fdata = pseudo.fdata,
                         pdata = pseudo.pdata)
  return(list(truep = true.prop, pseudo_eset = pseudo_eset))
}

#' Create an ExpressionSet object
#'
#' This function creates an ExpressionSet object which is a container for
#' storing gene expression data along with related experimental data.
#' An ExpressionSet object includes an expression matrix, feature data (fdata)
#' and phenotype data (pdata).
#'
#' @param exprs A matrix or data frame of expression values. Rows correspond to
#'  features (e.g., genes) and columns correspond to samples.
#' @param fdata A data frame or matrix of feature data. Each row corresponds to
#'  feature in the expression set and columns correspond to feature variables or
#'  a annotations. The row names should match the row names of the `exprs`
#'  parameter.
#' @param pdata A data frame or matrix of phenotype data. Each row corresponds
#' to a sample in the expression set and columns correspond to phenotype
#' variables. The row names should match the column names of the `exprs`
#' parameter.
#'
#' @return An ExpressionSet object which includes the expression data (`exprs`),
#' phenotype data (`pdata`), and feature data (`fdata`).
#'
#' @import Biobase
#' @import Matrix
#'
#' @export
getESET <- function (exprs, fdata, pdata)
{
  pdata <- as.data.frame(pdata)
  fdata <- as.data.frame(fdata)
  exprs <- Matrix::Matrix(exprs, sparse = TRUE)
  rownames(pdata) <- colnames(exprs)
  rownames(fdata) <- rownames(exprs)
  eset <- Biobase::ExpressionSet(exprs, Biobase::AnnotatedDataFrame(pdata),
                                 Biobase::AnnotatedDataFrame(fdata))
}
