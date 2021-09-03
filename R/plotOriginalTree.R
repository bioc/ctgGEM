#' Display Original ctgGEM Plots
#'
#' Displays the original plots created by the ctgGEM package and
#' stored in the [originalTrees] slot of a ctgGEMset object.
#'
#' @note In order to reproduce original plots, the
#'     respective package(s) must be installed.
#'
#' @param dataSet a ctgGEMset object
#' @param treeType the type of tree to display.  Must be one of
#'     \code{names(originalTrees(dataSet))}
#'
#' @return a ggplot2::ggplot object.
#' @export
#' @include ctgGEMset-class.R
#' @include ctgGEMset-methods.R
#' @importFrom graphics plot.new
#' @examples
#' # load HSMMSingleCell package
#' library(HSMMSingleCell)
#'
#' # load the data for TSCAN and monocle:
#' data(HSMM_expr_matrix)
#' data(HSMM_sample_sheet)
#' data(HSMM_gene_annotation)
#'
#' # construct a ctgGEMset
#' dataSet <- ctgGEMset(exprsData = HSMM_expr_matrix,
#'                         phenoData = HSMM_sample_sheet,
#'                         featureData = HSMM_gene_annotation)
#'
#' TSCANinfo(dataSet) <- "ENSG00000000003.10"
#'
#' # run generate_tree()
#' dataSet <- generate_tree(dataSet = dataSet, treeType = "TSCAN")
#'
#' # view names of original trees
#' names(originalTrees(dataSet))
#'
#' # plot original trees
#' plotOriginalTree(dataSet, "TSCANclustering")
#' plotOriginalTree(dataSet, "TSCANsingleGene")
plotOriginalTree <- function(dataSet, treeType) {
    stopifnot(is(dataSet, "ctgGEMset"))
    stopifnot(treeType %in% names(originalTrees(dataSet)))
    treeData <- originalTrees(dataSet)[[treeType]]

    # ANY CHANGES MADE IN THE FOLLOWING CODE MUST BE CHECKED FOR
    # COMPATIBILITY WITH THEIR RESPECTIVE make METHODS
# 
#     if (treeType == "cellTreeTopics") {
#         if (!requireNamespace("cellTree", quietly = TRUE)) {
#             stop(
#             "Package 'cellTree' is required to view an original cellTree plot.
#             See ctgGEM vignette for details on installing 'cellTree'",
#             call. = FALSE)
#         }
#         cellTree::ct.plot.topics(treeData)
#     } else if (treeType == "cellTreeGrouping") {
#         if (!requireNamespace("cellTree", quietly = TRUE)) {
#             stop(
#             "Package 'cellTree' is required to view an original cellTree plot.
#             See ctgGEM vignette for details on installing 'cellTree'",
#             call. = FALSE)
#         }
#         cellTree::ct.plot.grouping(treeData)
#     } else 
    if (grepl("sincell", treeType)) {
        igraph::plot.igraph(
            treeData,
            main = treeData$main,
            vertex.label = treeData$vertex.label,
            vertex.size = treeData$vertex.size,
            edge.color = treeData$edge.color,
            edge.width = treeData$edge.width,
            vertex.label.cex = treeData$vertex.label.cex,
            layout = treeData$layout
        )
    } else {
        # if it's not a special case, it's in ggplot format
        graphics::plot(treeData)
    }
}
