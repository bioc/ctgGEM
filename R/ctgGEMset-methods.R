#' @title Methods for the ctgGEMset class
#' @name ctgGEMset-methods
#' @aliases ctgGEMset-methods-methods
#'
#' @description  These methods operate on ctgGEMset objects. Please note that
#'    treeList<- and originalTrees<- are not intended to be
#'    called directly.
#' @param cs A ctgGEMset object
#' @param value
#'    \code{cellTreeInfo(cs)<-}: A character vector of the name of the
#'       column of \code{phenoData()} to use for grouping data
#'
#'    \code{monocleInfo(cs)<-}: The value to use as a named parameter for
#'       \code{generate_tree(treeType = "monocle")}
#'
#'    \code{TSCANinfo(cs)<-}: A character vector of the row name of a
#'       single gene in \code{exprsData()} to use for a single gene vs.
#'       pseudotime plot for \code{generate_tree(treeType = "TSCAN")}
#'
#'    \code{sincellInfo(cs)<-}: The value to use as a named parameter for
#'       sincell, used by \code{generate_tree(treeType = "sincell")}
#'
#' @param tt The type of tree being stored
#' @param pt The name of the \pkg{monocle} or \pkg{sincell} parameter to store
#' @return An updated ctgGEMset object, or the contents of a slot of the
#'      ctgGEMset object
#' @importFrom methods as is new
#' @examples
#' # load HSMMSingleCell package
#' library(HSMMSingleCell)
#'
#' # load the data
#' data(HSMM_expr_matrix)
#' data(HSMM_sample_sheet)
#' data(HSMM_gene_annotation)
#'
#' # construct a ctgGEMset
#' dataSet <- ctgGEMset(exprsData = HSMM_expr_matrix,
#'                         phenoData = HSMM_sample_sheet,
#'                         featureData = HSMM_gene_annotation)
#'
#' cellTreeInfo(dataSet) <- "Hours"
#' cellTreeInfo(dataSet)
#'
#' monocleInfo(dataSet, "gene_id") <- "gene_short_name"
#' monocleInfo(dataSet, "cell_id_1") <- "MYF5"
#' monocleInfo(dataSet, "cell_id_2") <- "ANPEP"
#' monocleInfo(dataSet, "ex_type") <- "FPKM"
#' monocleInfo(dataSet)
#'
#' TSCANinfo(dataSet) <- "ENSG00000000003.10"
#' TSCANinfo(dataSet)
#'
#' sincellInfo(dataSet, "method") <- "classical-MDS"
#' sincellInfo(dataSet, "MDS.distance") <- "spearman"
#' sincellInfo(dataSet, "clust.method") <- "k-medoids"
#' sincellInfo(dataSet)
#'
#' # The following two examples will return empty lists, since no trees
#' # have been generated on this ctgGEMset
#' trees <- treeList(dataSet)
#' originalTrees <- originalTrees(dataSet)
NULL

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export
cellTreeInfo <- function(cs) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@cellTreeInfo
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export
`cellTreeInfo<-` <- function(cs, value) {
    stopifnot(is(cs, "ctgGEMset"))
    if (length(value) != 1 || !is(value, "character") ||
        !value %in% colnames(colData(cs))) {
        stop(
            "cellTreeInfo must be a character vector of length 1,",
            " and a column name in phenoData()."
        )
    }
    cs@cellTreeInfo <- value
    cs
    }

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export
monocleInfo <- function(cs) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@monocleInfo
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export
`monocleInfo<-` <- function(cs, pt, value) {
    stopifnot(is(cs, "ctgGEMset"))
    valid_pt_types <- c("gene_id", "cell_id_1", "cell_id_2", "ex_type")
    valid_ex_type <- c("UMI", "TC", "FPKM", "TPM", "LTFPKM", "LTTPM")
    if(!(pt %in% valid_pt_types)){
        stop(pt, " is not a valid",
             " parameter type for monocleInfo.")
    }
    if(pt == "gene_id"){
        if(!(value %in% colnames(rowData(cs)))){
            stop("gene_id must be a column name in featureData().")
        }
    }
    if(pt == "cell_id_1" || pt == "cell_id_2"){
        if(is.null(monocleInfo(cs)[["gene_id"]])){
            stop("gene_id must be set before setting cell_id.")
        }
        valid_cells <- as.character(rowData(cs)[[monocleInfo(cs)[["gene_id"]]]])
        if(!(value %in% valid_cells)){
            stop("cell_id must be a gene ID found in the column set as gene_id")
        }
    }
    if(pt == "ex_type" && !(value %in% valid_ex_type)){
        stop("ex_type must be one of 'UMI', 'TC', 'FPKM', 'TPM', 'LTFPKM',",
                "'LTTPM'")
    }
    cs@monocleInfo[[pt]] <- value
    cs
    }


#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export

TSCANinfo <- function(cs) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@TSCANinfo
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export
`TSCANinfo<-` <- function(cs, value) {
    stopifnot(is(cs, "ctgGEMset"))
    if (length(value) != 1 || !is(value, "character") ||
        !(value %in% row.names(assay(cs)))) {
        stop("TSCANinfo must be a character vector of length 1, and a row name
                in the gene expression matrix.")
    }
    cs@TSCANinfo <- value
    cs
    }

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export

sincellInfo <- function(cs) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@sincellInfo
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export

`sincellInfo<-` <- function(cs, pt, value) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@sincellInfo[[pt]] <- value
    cs
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export

treeList <- function(cs) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@treeList
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @keywords internal

`treeList<-` <- function(cs, tt, value) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@treeList[[tt]] <- value
    cs
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @export

originalTrees <- function(cs) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@originalTrees
}

#' @rdname ctgGEMset-methods
#' @aliases ctgGEMset-methods
#' @keywords internal

`originalTrees<-` <- function(cs, tt, value) {
    stopifnot(is(cs, "ctgGEMset"))
    cs@originalTrees[[tt]] <- value
    cs
}
