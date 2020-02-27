#' A Cell Tree Generating Function using sincell
#'
#' This function, called by \code{\link{generate_tree}}, creates visualizations
#' using the data using the \pkg{sincell} package. This function utilizes code
#' from the sincell vignette R script.
#'
#' @param dataSet a ctgGEMset object
#' @param outputDir the directory where output should be saved, defaults to
#' the temporary location returned by \code{tempdir()}
#' @return an updated ctgGEMset object
#' @keywords internal
#' @import Biobase
#' @importFrom utils head tail
#' @importFrom methods is

makeSincell <- function(dataSet,outputDir = tempdir()){
    if (!requireNamespace("sincell", quietly = TRUE)) {
        stop(
            "Package 'sincell' is required for treeType = 'sincell'",
            "but is not installed.  See vignette for details on installing",
            "'sincell'",
            call. = FALSE
        )
    }
    # SO <- sincell::sc_InitializingSincellObject(exprs(dataSet))
    SO <- sincell::sc_InitializingSincellObject(SummarizedExperiment::assay(dataSet))
    params <- sincellInfo(dataSet)
    #get cell distances
    oned_dist_methods <- c("euclidean", "cosine", "pearson", "spearman",
                            "L1", "MI")
    md_dist_methods <- c("PCA", "ICA", "tSNE", "classical-MDS", "nonmetric-MDS")
    clust_methods <- c("max-distance", "percent", "knn", "k-medoids", "ward.D",
                        "ward.D2", "single", "complete", "average", "mcquitty",
                        "median", "centroid")
    if (!is.null(params[["method"]]) &&
        params[["method"]] %in% md_dist_methods) {
        # cosine distance not supported for dimension reduction
        if (!is.null(params[["MDS.distance"]])) {
            SO <- sincell::sc_DimensionalityReductionObj(
                SO, method = params[["method"]],
                MDS.distance = params[["MDS.distance"]])
        } else {
            SO <- sincell::sc_DimensionalityReductionObj(
                    SO, method = params[["method"]])
        }
    } else if (!is.null(params[["method"]])
                && params[["method"]] %in% oned_dist_methods) {
        SO <- sincell::sc_distanceObj(SO, method = params[["method"]])
    } else {
        message("'method' in sincellInfo missing or invalid.")
        message("Computing dimensionality reduction and distance using defaults.")
        message("See ctgGEM vignette for more information.")
        SO <- sincell::sc_DimensionalityReductionObj(SO)
    }
    #clustering
    if (!is.null(params[["clust.method"]])
        && params[["clust.method"]] %in% clust_methods) {
        SO <-
            sincell::sc_clusterObj(SO, clust.method = params[["clust.method"]])
    } else {
        message("'clust.method' in sincellInfo missing or invalid.")
        message("Computing clustering using default.")
        message("See ctgGEM vignette for more information.")
        SO <- sincell::sc_clusterObj(SO)
    }
    #build the trees
    # Minimum Spanning Tree (MST)
    SO <- sincell::sc_GraphBuilderObj(SO, graph.algorithm = "MST",
                                        graph.using.cells.clustering = TRUE)
    cellstateHierarchy_MST <- SO[["cellstateHierarchy"]]
    # Maximum Similarity Spanning Tree (SST)
    SO <- sincell::sc_GraphBuilderObj(SO, graph.algorithm = "SST",
                                        graph.using.cells.clustering = TRUE)
    cellstateHierarchy_SST <- SO[["cellstateHierarchy"]]
    # Iterative Mutual Clustering Graph (IMC)
    SO <- sincell::sc_GraphBuilderObj(SO, graph.algorithm = "IMC")
    cellstateHierarchy_IMC <- SO[["cellstateHierarchy"]]

    # set plotting parameters
    igraph::V(cellstateHierarchy_MST)$vertex.size <- 5
    igraph::V(cellstateHierarchy_SST)$vertex.size <- 5
    igraph::V(cellstateHierarchy_IMC)$vertex.size <- 5

    igraph::V(cellstateHierarchy_MST)$vertex.label.cex <- 0.2
    igraph::V(cellstateHierarchy_SST)$vertex.label.cex <- 0.2
    igraph::V(cellstateHierarchy_IMC)$vertex.label.cex <- 0.2

    igraph::V(cellstateHierarchy_MST)$vertex.label <-
        colnames(SO[["expressionmatrix"]])
    igraph::V(cellstateHierarchy_SST)$vertex.label <-
        colnames(SO[["expressionmatrix"]])
    igraph::V(cellstateHierarchy_IMC)$vertex.label <-
        colnames(SO[["expressionmatrix"]])

    igraph::E(cellstateHierarchy_MST)$edge.color <- "black"
    igraph::E(cellstateHierarchy_SST)$edge.color <- "black"
    igraph::E(cellstateHierarchy_IMC)$edge.color <- "black"

    igraph::E(cellstateHierarchy_MST)$edge.width <- 2
    igraph::E(cellstateHierarchy_SST)$edge.width <- 2
    igraph::E(cellstateHierarchy_IMC)$edge.width <- 2

    cellstateHierarchy_MST$layout <- igraph::layout.kamada.kawai
    cellstateHierarchy_SST$layout <- igraph::layout.kamada.kawai
    cellstateHierarchy_IMC$layout <- igraph::layout.kamada.kawai

    cellstateHierarchy_MST$main <- "MST"
    cellstateHierarchy_SST$main <- "SST"
    cellstateHierarchy_IMC$main <- "IMC"

    #layout.graph=igraph::layout.kamada.kawai;
    # par(bty="o",xaxs="i",yaxs="i",cex.axis=mycex-0.2,cex.main=mycex,cex.
    #     lab=mycex,las=1,mar=c(5.3,5.3,2.9,1.6),oma=c(1,1,2,10))
    #format the filename
    filename <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    
    # plot and save sincell MST
    if(is.null(outputDir)){
        fn <- tempfile(paste0(filename,"_sincellMST"),
                       tmpdir = file.path(tempdir(),"CTG-Output","Plots"),fileext=".png")
        grDevices::png(filename = fn)
    } else {
        #open png writer
        grDevices::png(filename = file.path(outputDir,"CTG-Output","Plots",
                                            paste0(filename,"_sincellMST.png")))
    }
    igraph::plot.igraph(cellstateHierarchy_MST,
        main = cellstateHierarchy_MST$main,
        vertex.label = cellstateHierarchy_MST$vertex.label,
        vertex.size = cellstateHierarchy_MST$vertex.size,
        edge.color = cellstateHierarchy_MST$edge.color,
        edge.width = cellstateHierarchy_MST$edge.width,
        vertex.label.cex = cellstateHierarchy_MST$vertex.label.cex,
        layout = cellstateHierarchy_MST$layout)
    grDevices::dev.off()
    
    # plot and save sincell SST
    if(is.null(outputDir)){
        fn <- tempfile(paste0(filename,"_sincellSST"),
                       tmpdir = file.path(tempdir(),"CTG-Output","Plots"),fileext=".png")
        grDevices::png(filename = fn)
    } else {
        #open png writer
        grDevices::png(filename = file.path(outputDir,"CTG-Output","Plots",
                                            paste0(filename,"_sincellSST.png")))
    }
    igraph::plot.igraph(cellstateHierarchy_SST,
        main = cellstateHierarchy_SST$main,
        vertex.label = cellstateHierarchy_SST$vertex.label,
        vertex.size = cellstateHierarchy_SST$vertex.size,
        edge.color = cellstateHierarchy_SST$edge.color,
        edge.width = cellstateHierarchy_SST$edge.width,
        vertex.label.cex = cellstateHierarchy_SST$vertex.label.cex,
        layout = cellstateHierarchy_SST$layout)
    grDevices::dev.off()
    
    # plot and save sincell IMC
    if(is.null(outputDir)){
        fn <- tempfile(paste0(filename,"_sincellIMC"),
                       tmpdir = file.path(tempdir(),"CTG-Output","Plots"),fileext=".png")
        grDevices::png(filename = fn)
    } else {
        #open png writer
        grDevices::png(filename = file.path(outputDir,"CTG-Output","Plots",
                                            paste0(filename,"_sincellIMC.png")))
    }
    igraph::plot.igraph(cellstateHierarchy_IMC,
        main = cellstateHierarchy_IMC,
        vertex.label = cellstateHierarchy_IMC$vertex.label,
        vertex.size = cellstateHierarchy_IMC$vertex.size,
        edge.color = cellstateHierarchy_IMC$edge.color,
        edge.width = cellstateHierarchy_IMC$edge.width,
        vertex.label.cex = cellstateHierarchy_IMC$vertex.label.cex,
        layout = cellstateHierarchy_IMC$layout)
    grDevices::dev.off()
    # ANY CHANGES MADE IN THE FOLLOWING 3 LINES OF CODE MUST BE CHECKED FOR
    # COMPATIBILITY WITH plotOriginalTree
    originalTrees(dataSet, "sincellMST") <- cellstateHierarchy_MST
    originalTrees(dataSet, "sincellSST") <- cellstateHierarchy_SST
    originalTrees(dataSet, "sincellIMC") <- cellstateHierarchy_IMC
    # convert data to standard cell tree format
    MST <- sincell2CTF(cellstateHierarchy_MST, filename, "MST", outputDir)
    SST <- sincell2CTF(cellstateHierarchy_SST, filename, "SST", outputDir)
    IMC <- sincell2CTF(cellstateHierarchy_IMC, filename, "IMC", outputDir)
    # store simplified igraph objects
    treeList(dataSet, "sincellMST") <- tree2igraph(MST)
    treeList(dataSet, "sincellSST") <- tree2igraph(SST)
    treeList(dataSet, "sincellIMC") <- tree2igraph(IMC)
    dataSet
}

# This helper function converts a sincell tree to the
# standard cell tree format and writes its SIF file.

sincell2CTF <- function(tree, filename, treeType, outputDir = NULL) {
    cellEdges <- igraph::ends(tree, es = igraph::E(tree))
    if (treeType == "MST") {
        relationshipType <- "minimum span"
    } else if (treeType == "SST") {
        relationshipType <- "maximum similarity"
    } else {
        relationshipType <- "iterative mutual clustering"
    }
    relationships <-
        paste0(cellEdges[, 1], '\t', relationshipType, '\t',
                cellEdges[, 2])
    # write these relationships to file
    if(is.null(outputDir)){
        fileName <- tempfile(paste0(filename,"_SIN_",treeType,"_CTF"),
                             tmpdir = file.path(tempdir(),"CTG-Output","SIFs"),fileext=".sif")
    } else {
        fileName <- file.path(outputDir,"CTG-Output","SIFs",
                              paste0(filename,"_SIN_",treeType,"_CTF.sif"))
    }
    # fullFileName <- file.path(outputDir,"CTG-Output","SIFs",
    #                           paste0(filename,"_SIN_",treeType,"_CTF.sif"))
    write(relationships, fileName)
    relationships
}
