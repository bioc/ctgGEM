#' A Cell Tree Generating Function using cellTree
#'
#' This function, called by \code{\link{generate_tree}}, creates a cell tree
#' using the data and the \pkg{cellTree} package. This function utilizes code
#' from the cellTree vignette R script.
#'
#' @param dataSet a ctgGEMset object
#' @return an updated ctgGEMset object
#' @keywords internal
#' @importFrom igraph vertex_attr ends E

makeCellTree <- function(dataSet) {
    if (!requireNamespace("cellTree", quietly = TRUE)) {
        stop(
            "Package 'cellTree' is required for treeType = 'cellTree',
            but is not installed.  See vignette for details on installing
            'cellTree'",
            call. = FALSE
        )
    }
    if (length(cellTreeInfo(dataSet)) == 0) {
        groupInfo <- NULL
    } else {
        # groupInfo <- as.character(
        #   Biobase::pData(dataSet)[[cellTreeInfo(dataSet)]])
        groupInfo <- as.character(
            SummarizedExperiment::colData(dataSet)[[cellTreeInfo(dataSet)]]
        )
    }
    
    # d <- as.matrix(Biobase::exprs(dataSet))
    d <- as.matrix(SummarizedExperiment::assay(dataSet))
    # compute lda results, using maptpx method and finding optimal number of
    # topics, between 2 and 15
    lda.results <- cellTree::compute.lda(d)
    # create backbone tree
    b.tree <-
        cellTree::compute.backbone.tree(lda.results, grouping = groupInfo)
    # format the filename
    filename <- as.character(Sys.time())
    filename <- gsub("/", "-", filename)
    filename <- gsub(":", "-", filename)
    filename <- gsub(" ", "_", filename)
    grDevices::png(filename = paste0("./CTG-Output/Plots/", filename,
                                     "_cellTreeTopics.png"))
    cellTree::ct.plot.topics(b.tree)
    grDevices::dev.off()
    # ANY CHANGES MADE IN THE FOLLOWING LINE OF CODE MUST BE CHECKED FOR
    # COMPATIBILITY WITH plotOriginalTree
    originalTrees(dataSet, "cellTreeTopics") <- b.tree
    if(!is.null(groupInfo)){
        grDevices::png(filename = paste0("./CTG-Output/Plots/", filename,
                                         "_cellTreeGrouping.png"))
        cellTree::ct.plot.grouping(b.tree)
        grDevices::dev.off()
        # ANY CHANGES MADE IN THE FOLLOWING LINE OF CODE MUST BE CHECKED FOR
        # COMPATIBILITY WITH plotOriginalTree
        originalTrees(dataSet, "cellTreeGrouping") <- b.tree
    }
    # convert data to standard cell tree format
    tree <- CT2CTF(b.tree, filename)
    treeList(dataSet, "cellTree") <- tree2igraph(tree)
    dataSet
}

# This helper function converts a cellTree cell tree to the standard cell
# tree format and writes its SIF file.

CT2CTF <- function(bTree, timeStamp) {
    # get the attributes of the backbone tree
    cellAttr <- vertex_attr(bTree)
    # get the edges (labeled with numbered ids) in the backbone tree
    cellEdges <- ends(bTree, es = igraph::E(bTree))
    # get the index of the edges that have a vertebrae to vertebrae connection
    v2v <- which(cellAttr$is.backbone[cellEdges[, 1]]
                    & cellAttr$is.backbone[cellEdges[, 2]])
    # relate the cells that form the backbone of the backbone tree
    relationships <- paste0(cellAttr$cell.name[cellEdges[v2v, 1]],
                            "\tvertebrae to vertebrae\t",
                            cellAttr$cell.name[cellEdges[v2v, 2]])
    relationships2 <- paste0(cellAttr$cell.name[cellEdges[, 1]],
                                "\tlow Chi-Square distance\t",
                                cellAttr$cell.name[cellEdges[, 2]])
    relationships <- append(relationships, relationships2)
    # write these relationships to file
    fileName <- paste0("./CTG-Output/SIFs/", timeStamp, "_CT_CTF.sif")
    write(relationships, fileName)
    relationships
}
