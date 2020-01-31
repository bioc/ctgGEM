#' Function to Generate Cell Trees
#'
#' This function builds a cell hierarchy tree of a chosen supported type with
#' a given data set, contained in a ctgGEMset object.  Different tree
#' types require data from corresponding slots of the ctgGEMset object.
#' See vignette for examples, usage details, and instructions on building a
#' ctgGEMset object.
#' @param dataSet the ctgGEMset object for creating the cell tree
#' @param treeType the type of tree generated
#' @return An updated ctgGEMset object.  The generated tree is placed in
#'     \code{@@treeList[treeType]} slot, and can be accessed via
#'     \code{treeList(dataSet)$treeType}.  The function also creates a
#'     directory named "treeType-Output" and writes the plot(s) of the
#'     generated tree(s) and its SIF file to that directory.
#'
#' @keywords cell tree
#' @export
#' @include ctgGEMset-class.R
#' @include ctgGEMset-methods.R
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
#' dataSet <- newctgGEMset(exprsData = HSMM_expr_matrix,
#'                         phenoData = HSMM_sample_sheet,
#'                         featureData = HSMM_gene_annotation)
#'
#' TSCANinfo(dataSet) <- "ENSG00000000003.10"
#'
#' # run generate_tree()
#' dataSet <- generate_tree(dataSet = dataSet, treeType = "TSCAN")

generate_tree <- function(dataSet, treeType) {
    stopifnot(is(dataSet, "ctgGEMset"))
    # get the matrix containing packages and corresponding functions
    method_matrix <- make_method_matrix()
    # find the correct package to use
    pack <- which(method_matrix[1, ] == treeType)
    # create a directory and subdirectories for output, if they doesn't exist
    if (!utils::file_test(op = "-d", "./CTG-Output")) {
        dir.create("./CTG-Output")
        dir.create("./CTG-Output/Plots")
        dir.create("./CTG-Output/SIFs")
    }
    # get the corresponding function
    func <- get(method_matrix[2, pack])
    # execute the function
    dataSet <- func(dataSet)
    dataSet
}

#' A Function That Creates A Predefined Matrix Of Usable Packages
#'
#' This function constructs the predefined matrix that contains all possible
#' cell tree building packages and their corresponding functions
#' @keywords internal
#' @return a matrix containing all possible packages and their corresponding
#'     function

make_method_matrix <- function() {
    # create a matrix of the packages and corresponding function calls
    method_matrix = matrix(
        c(
            'cellTree',
            'TSCAN',
            'monocle',
            'destiny',
            'sincell',
            'makeCellTree',
            'makeTSCAN',
            'makeMonocle',
            'makeDestiny',
            'makeSincell'
        ),
        nrow = 2,
        ncol = 5,
        byrow = TRUE,
        dimnames = list(c('package', 'function'))
    )
    # return the method matrix to the caller
    return(method_matrix)
}
