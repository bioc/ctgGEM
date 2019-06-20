#' Convert Cell Relationships in a SIF object to igraph Format
#'
#' This function converts cell relationships in a cell tree into igraph format
#' @param relationships character vector of a cell tree in SIF format
#' @return a cell tree in igraph format
#' @importFrom igraph make_empty_graph add_edges set_edge_attr set_vertex_attr
#' @importFrom igraph %>%
#' @keywords internal

tree2igraph <- function(relationships){
    # turn the lines into a matrix of relationships with columns 1 and 3 being
    # related items and column 2 being the relationship type
    relationshipsMatrix <- matrix(unlist(strsplit(relationships, "\t")),
                                    ncol=3, byrow = TRUE)
    # extract the relationship types from the matrix
    relationshipTypes <- relationshipsMatrix[,-c(1,3)]
    # extract the related items into a vector by row
    relatedItems <- as.vector(t(relationshipsMatrix[,-2]))
    # get the unique items from related items to get the indices
    vertexIDs <- unique(relatedItems)
    # identify the related items by number instead of name
    relatedIDs <- vector(mode = "integer", length = length(relatedItems))
    relatedIDs <- vapply(relatedItems, function(x){
        return(which(vertexIDs %in% x))
    }, 1L)
    cellIgraph <- make_empty_graph(n = length(vertexIDs)) %>%
                    add_edges(relatedIDs)  %>%
                    set_edge_attr("label", value = relationshipTypes) %>%
                    set_vertex_attr("name", value = vertexIDs)
    cellIgraph
}
