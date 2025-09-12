## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", width = 80)
knitr::opts_knit$set(global.par = TRUE)

## ----echo = FALSE-------------------------------------------------------------
set.seed(5)
oldmar <- par()$mar
par(mar = c(0, 0, 0, 0) + 0.1)

## -----------------------------------------------------------------------------
library(backbone)

## -----------------------------------------------------------------------------
library(igraph)

## -----------------------------------------------------------------------------
W <- matrix(c(0,10,10,10,10,75,0,0,0,0,  #Adjacency matrix of example network
              10,0,1,1,1,0,0,0,0,0,
              10,1,0,1,1,0,0,0,0,0,
              10,1,1,0,1,0,0,0,0,0,
              10,1,1,1,0,0,0,0,0,0,
              75,0,0,0,0,0,100,100,100,100,
              0,0,0,0,0,100,0,10,10,10,
              0,0,0,0,0,100,10,0,10,10,
              0,0,0,0,0,100,10,10,0,10,
              0,0,0,0,0,100,10,10,10,0),10)
W <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)  #Convert to an igraph object

## -----------------------------------------------------------------------------
plot(W, edge.width = sqrt(E(W)$weight), vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- backbone_from_weighted(W, model = "global", parameter = 0)
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- backbone_from_weighted(W, model = "global", parameter =  mean(E(W)$weight))
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- backbone_from_weighted(W, model = "disparity", alpha = 0.05)
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),  #Incidence matrix of example bipartite network
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
B <- graph_from_biadjacency_matrix(B)  #Convert to an igraph object
plot(B, layout = layout_as_bipartite(B), vertex.label = NA, vertex.size = 1)

## -----------------------------------------------------------------------------
P <- bipartite_projection(B, which = "true")
plot(P, vertex.label = NA, edge.width = sqrt(E(P)$weight), edge.color = rgb(0,0,0,.1))

## -----------------------------------------------------------------------------
bb <- backbone_from_projection(B, model = "sdsm", alpha = 0.05)
plot(bb, vertex.label = NA)

## ----fig.width=8, out.width="75%"---------------------------------------------
par(mfrow = c(1, 2), mar = c(0,0,2,0))  #Display plots side-by-side
U <- sample_sbm(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))  #Example unweighted network
plot(U, main = "Original Network", vertex.label = NA)
bb <- backbone_from_unweighted(U, model = "lspar", parameter = 0.5)
plot(bb, main = "L-Spar Backbone", vertex.label = NA)  #Communities are clearly visible

## ----fig.width=8, out.width="75%"---------------------------------------------
par(mfrow = c(1, 2), mar = c(0,0,2,0))  #Display plots side-by-side
U <- sample_pa(n = 60, m = 3, directed = FALSE)  #Example unweighted network
V(U)$name <- ""
V(U)$name[which(degree(U)>=sort(degree(U), decreasing = TRUE)[5])] <- LETTERS[1:5]  #Label five highest-degree nodes
plot(U, vertex.size = degree(U), vertex.label.family = "sans", vertex.label.font = 2, main = "Original Network")
bb <- backbone_from_unweighted(U, model = "degree", parameter = 0.25)
plot(bb, vertex.size = degree(bb), vertex.label.family = "sans", vertex.label.font = 2, main = "Local Degree Backbone")

## -----------------------------------------------------------------------------
B <- matrix(sample(c(0,1), 18, replace = TRUE), 3, 6)  #A simple incidence matrix
bb <- backbone_from_projection(B, model = "sdsm", alpha = 0.05, return = "everything")  #Extract backbone, return everything
bb$bipartite  #The original bipartite network
bb$projection  #The weighted projection
bb$pvalues  #The edgewise p-values
bb$backbone  #The backbone (as a matrix)
bb$narrative  #Narrative description
bb$call  #Function call

## ----fig.width=8, out.width="75%"---------------------------------------------
par(mfrow = c(1, 2), mar = c(0,0,2,0))  #Display plots side-by-side
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),  #Incidence matrix of example bipartite network
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
B <- graph_from_biadjacency_matrix(B)  #Convert to an igraph object
bb1 <- backbone_from_projection(B, model = "sdsm", alpha = 0.05)  #Extract non-signed backbone
layout <- layout_nicely(bb1)
plot(bb1, vertex.labels = NA, layout = layout, main = "Non-signed")
bb2 <- backbone_from_projection(B, model = "sdsm", alpha = 0.1, signed = TRUE)  #Extract signed backbone
E(bb2)$color <- "green"
E(bb2)$color[which(E(bb2)$sign==-1)] <- "red"
plot(bb2, vertex.labels = NA, layout = layout, main = "Signed")

## ----fig.width=8, out.width="75%"---------------------------------------------
par(mfrow = c(1, 2), mar = c(0,0,2,0))  #Display plots side-by-side
U <- sample_gnp(100, .25)
plot(U, vertex.label = NA, main = "Original Network")
bb <- backbone_from_unweighted(U, model = "custom", escore = "jaccard", normalize = "none", filter = "proportion", parameter = 0, umst = TRUE)
plot(bb, vertex.label = NA, main = "Custom Backbone")

## -----------------------------------------------------------------------------
mat <- rbind(c(1,0,0), c(0,0,1), c(0,1,1))
mat
fastball(mat)

## -----------------------------------------------------------------------------
mat <- rbind(c(1,0,0), c(0,0,1), c(0,1,1))
bicm(mat)

