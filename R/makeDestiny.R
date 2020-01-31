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
        assayData = exprs(dataSet),
        phenoData = AnnotatedDataFrame(pData(dataSet)),
        featureData = AnnotatedDataFrame(fData(dataSet))
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

    # convert data to standard cell tree format
    tree <- destiny2CTF(dpt, filename)
    treeList(dataSet, "destiny") <- tree2igraph(tree)
    dataSet
}

# The following function was 'lifted' from
# https://github.com/sirusb/destiny/blob/master/R/dpt-plotting.r
# It is used by the destiny package (and here) to compute pseudotime ordering,
# but it is not exported by the destiny package

dpt_for_branch <- function(dpt, branch_id) {
    branch_idx <- dpt@branch[, 1L] == branch_id
    stopifnot(any(branch_idx))
    tip_cells <- which(branch_idx & dpt@tips[, 1L])
    if (length(tip_cells) == 0L) tip_cells <- which(branch_idx)
    dpt[tip_cells[[1L]], ]
}


# This helper function converts a destiny Diffusion PseudoTime tree to the
# standard cell tree format and writes its SIF file.

destiny2CTF <- function(dpt, filename) {
    # reconstruct diffusion pseudotime from destiny's plotting methods
    root <- min(dpt@branch[, 1], na.rm = TRUE)
    pt_vec <- dpt_for_branch(dpt, root) # vector containing pseudotime
    names(pt_vec) <- rownames(destiny::as.data.frame(dpt))
    sorted_pt_vec <- sort(pt_vec)
    cellOrder <- names(sorted_pt_vec)
    nextCellOrder <- c(tail(cellOrder,-1), NA)
    relationships <-
        paste0(cellOrder, "\tpsuedotime\t", nextCellOrder)
    # remove the last element because it contains no second cell
    relationships <- head(relationships,-1)
    # write these relationships to file
    write(relationships,
        paste0("./CTG-Output/SIFs/", filename, "_DPT_CTF.sif"))
    relationships
}


