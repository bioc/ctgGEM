#' A Cell Tree Generating Function using destiny
#'
#' This function, called by \code{\link{generate_tree}}, creates visualizations
#' using the data using the \pkg{destiny} package. This function utilizes code
#' from the destiny and dpt vignette R scripts.
#'
#' @param dataSet a ctgGEMset object
#' @return an updated ctgGEMset object
#' @keywords internal
#' @import Biobase
#' @importFrom utils head tail
#' @importFrom methods is
#' @importFrom grDevices palette
#' @importFrom graphics title
#'
makeDestiny <- function(dataSet) {
    if (!requireNamespace("destiny", quietly = TRUE)) {
        stop(
            "Package 'destiny' is required for treeType = 'destiny',
            but is not installed.  See vignette for details on installing
            'destiny'",
            call. = FALSE
        )
    }
    es <- ExpressionSet(
        assayData = assayData(dataSet),
        phenoData = phenoData(dataSet),
        featureData = featureData(dataSet)
    )
    dpt <- destiny::DiffusionMap(es)
    dpt <- destiny::DPT(dpt)
    #format the filename
    filename <- as.character(Sys.time())
    filename <- gsub("/", "-", filename)
    filename <- gsub(":", "-", filename)
    filename <- gsub(" ", "_", filename)
    # configure color palette
    palette(destiny::cube_helix(6))
    #open png writer
    grDevices::png(filename = paste0(
        "./CTG-Output/Plots/",
        filename,
        "_destinyDiffusionMap.png"
    ))
    # generate the plot for diffusionmap
    destiny::plot(dpt@dm, main = "DiffusionMap")
    #close the writing device
    grDevices::dev.off()

    # generate the plot for diffusionmap
    destiny::plot(dpt, main = "DPT")
    ggplot2::ggsave(filename = paste0("./CTG-Output/Plots/", filename,
                                        "_destinyDPT.png"))
    # store the original plots, using a placeholder for DM since it's in DPT
    # ANY CHANGES MADE IN THE FOLLOWING 2 LINES OF CODE MUST BE CHECKED FOR
    # COMPATIBILITY WITH plotOriginalTree
    originalTrees(dataSet, "destinyDM") <- 1
    originalTrees(dataSet, "destinyDPT") <- dpt

    # convert data to standard cell tree format (Not ready yet)
    # tree <- destiny2CTF(dpt, filename)
    # treeList(dataSet, "destiny") <- tree2igraph(tree)
    dataSet
}

# This helper function converts a destiny Diffusion PseudoTime tree to the
# standard cell tree format and writes its SIF file. Currently a useless stub.

destiny2CTF <- function(tree, filename) {
    NULL
}
