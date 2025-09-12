## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", width = 80)
knitr::opts_knit$set(global.par = TRUE)

## ----echo = FALSE-------------------------------------------------------------
set.seed(5)
oldmar <- par()$mar
par(mar = c(0, 0, 2, 0) + 0.1)

## -----------------------------------------------------------------------------
library(backbone)
library(igraph)

## -----------------------------------------------------------------------------
data(senate108)
senate108

## -----------------------------------------------------------------------------
as_biadjacency_matrix(senate108)[1:5,1:5]

## -----------------------------------------------------------------------------
rowSums(as_biadjacency_matrix(senate108))[1:5]
colSums(as_biadjacency_matrix(senate108))[1:5]

## -----------------------------------------------------------------------------
projection <- bipartite_projection(senate108, which = "false")
V(projection)$name <- V(projection)$last #Use only last names as node labels
as_adjacency_matrix(projection, attr = "weight")[1:5,1:5]
plot(projection, vertex.label = NA, vertex.frame.color = NA, vertex.size = 3, edge.width = E(projection)$weight^.1, edge.color = rgb(0,0,0,.1), main = "Weighted Projection")

## -----------------------------------------------------------------------------
bb1 <- backbone_from_projection(senate108, model = "sdsm", alpha = 0.05, narrative = TRUE)
layout <- layout_nicely(bb1)  #Get layout for backbone (we'll use it later)
plot(bb1, vertex.label = NA, vertex.frame.color = NA, vertex.size = 3, edge.color = rgb(0,0,0,.1), layout = layout, main = "SDSM Backbone")

## -----------------------------------------------------------------------------
bb1_signed <- backbone_from_projection(senate108, model = "sdsm", alpha = 0.1, narrative = TRUE, signed = TRUE)
E(bb1_signed)$color <- "green"  #Make all edges green
E(bb1_signed)$color[which(E(bb1_signed)$sign == -1)] <- rgb(1,0,0,.05)  #Make negative edges transparent red
plot(bb1_signed, vertex.label = NA, vertex.frame.color = NA, vertex.size = 3, layout = layout, main = "SDSM Signed Backbone")

## -----------------------------------------------------------------------------
bb2 <- backbone_from_weighted(projection, model = "disparity", alpha = 0.2, narrative = TRUE)
plot(bb2, vertex.label = NA, vertex.frame.color = NA, vertex.size = 3, edge.color = rgb(0,0,0,.1), main = "Disparity Filter Backbone")

## -----------------------------------------------------------------------------
unweighted <- delete_edges(projection, which(E(projection)$weight < 25))  #Delete low-weight edges
unweighted <- delete_edge_attr(unweighted, "weight")  #Delete edge weights to obtain an unweighted network
unweighted <- delete_vertices(unweighted, which(degree(unweighted) < 1))  #Delete isolated nodes
edge_density(unweighted)  #Compute density
plot(unweighted, vertex.label = NA, vertex.frame.color = NA, vertex.size = 3, edge.color = rgb(0,0,0,.25), main = "Unweighted Network")

## -----------------------------------------------------------------------------
bb3 <- backbone_from_unweighted(unweighted, model = "lspar", parameter = .5, narrative = TRUE)
plot(bb3, vertex.label = NA, vertex.frame.color = NA, vertex.size = 3, edge.color = rgb(0,0,0,.25), main = "L-Spar Backbone")

