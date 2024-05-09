## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
knitr::opts_knit$set(global.par = TRUE)

## ----echo = FALSE-------------------------------------------------------------
set.seed(5)
oldmar <- par()$mar
par(mar = c(0, 0, 0, 0) + 0.1)

## ----setup--------------------------------------------------------------------
library(backbone)

## -----------------------------------------------------------------------------
dat <- matrix(runif(100),10,10)  #Some data

backbone.suggest(dat)  #What should I do?

backbone <- backbone.suggest(dat, s = 0.05)  #Or, just do it

## -----------------------------------------------------------------------------
W <- matrix(c(0,10,10,10,10,75,0,0,0,0,
              10,0,1,1,1,0,0,0,0,0,
              10,1,0,1,1,0,0,0,0,0,
              10,1,1,0,1,0,0,0,0,0,
              10,1,1,1,0,0,0,0,0,0,
              75,0,0,0,0,0,100,100,100,100,
              0,0,0,0,0,100,0,10,10,10,
              0,0,0,0,0,100,10,0,10,10,
              0,0,0,0,0,100,10,10,0,10,
              0,0,0,0,0,100,10,10,10,0),10)

## -----------------------------------------------------------------------------
weighted <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)
plot(weighted, edge.width = sqrt(igraph::E(weighted)$weight), vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- global(W, upper = 0, class = "igraph")
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- global(W, upper = function(x)mean(x), class = "igraph")
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- disparity(W, alpha = 0.05, narrative = TRUE, class = "igraph")
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))

## -----------------------------------------------------------------------------
B[1:5,1:5]

## -----------------------------------------------------------------------------
rowSums(B)
colSums(B)

## -----------------------------------------------------------------------------
P <- B%*%t(B)
plot(igraph::graph_from_adjacency_matrix(P, mode = "undirected", diag = FALSE, weighted = TRUE), vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- sdsm(B, alpha = 0.075, narrative = TRUE, class = "igraph")

## -----------------------------------------------------------------------------
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
U.with.communities <- igraph::sbm.game(60, matrix(c(.75,.25,.25,.25,.75,.25,.25,.25,.75),3,3), c(20,20,20))
plot(U.with.communities, vertex.label = NA)

## -----------------------------------------------------------------------------
bb <- sparsify.with.lspar(U.with.communities, s = 0.6, narrative = TRUE)
plot(bb, vertex.label = NA)

## -----------------------------------------------------------------------------
U.with.hubs <- igraph::as.undirected(igraph::sample_pa(60, m = 3), mode = "collapse")
plot(U.with.hubs, vertex.size = igraph::degree(bb), vertex.label = NA) #A hairball

## -----------------------------------------------------------------------------
bb <- sparsify.with.localdegree(U.with.hubs, s = 0.3, narrative = TRUE)
plot(bb, vertex.size = igraph::degree(bb), vertex.label = NA)

## ----echo = TRUE, results = 'hide', warning = FALSE---------------------------
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
bb.object <- fdsm(B, alpha = NULL, trials = 1000)  #Backbone object containing edgewise p-values

## -----------------------------------------------------------------------------
bb1 <- backbone.extract(bb.object, alpha = 0.5, class = "igraph")  #Backbone extracted at alpha = 0.5
plot(bb1, vertex.label = NA)

## -----------------------------------------------------------------------------
bb2 <- backbone.extract(bb.object, alpha = 0.05, class = "igraph")  #Backbone extracted at alpha = 0.05
plot(bb2, vertex.label = NA)

## -----------------------------------------------------------------------------
mat <- rbind(c(1,0,0), c(0,0,1), c(0,1,1))
mat
fastball(mat)

## -----------------------------------------------------------------------------
mat <- rbind(c(1,0,0), c(0,0,1), c(0,1,1))
bicm(mat)

## ----echo = FALSE-------------------------------------------------------------
par(mar = oldmar)

