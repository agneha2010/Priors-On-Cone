library(rnhanesdata)
require(stats)
require(graphics)
library(splines2)
library(coneproj)
library(rcdd)
library(truncnorm)
library(truncdist)
library(mvtnorm)
library(R.utils)
library(Matrix)
library(igraph)
library(MCMCpack)
library(fda)
library(nnls)

###############model fitting

source("/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava_v1/fittingCompoSpikeFinal - Mixed effect Data.R")
source('/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava_v1/fittingCompo - Mixed effect Data.R')

Ti=24
t=1:24
i1 = 9
i2 = 20
A = matrix(0,nrow=Ti+2,ncol=Ti)
for(j in 1:(i1-2)){
  A[j,j] = t[(j+2)]-t[(j+1)]
  A[j,(j+1)] = t[j]-t[(j+2)]
  A[j,(j+2)] = t[j+1]-t[j]
}

for(j in (i2):(Ti-2)){
  A[j,j] = t[(j+2)]-t[(j+1)]
  A[j,(j+1)] = t[j]-t[(j+2)]
  A[j,(j+2)] = t[j+1]-t[j]
}

for(j in (i1-1):(i2-1)){
  A[j,j] = t[(j+1)]-t[(j+2)]
  A[j,(j+1)] = -t[j]+t[(j+2)]
  A[j,(j+2)] = -t[j+1]+t[j]
}
## added by me
#A[i2-1,(i2-2):i2] = - A[i2-1,(i2-2):i2]

A[(Ti-1),1] = -1.0
A[(Ti-1),2] = 1.0
A[Ti,1] = 1.0
A[(Ti+1),(Ti-1)] = 1.0
A[(Ti+1),Ti] = -1.0
A[(Ti+2),Ti] = 1.0

J    <- 24
knot <- J-3

BS1 <- bsplineS((1:Ti)/Ti, breaks = seq(0, 1, 1/knot))

A <- A %*% BS1

m = nrow(A)
### makeH uses Ax <= 0; so multiply by -1
B = -A
qux <- makeH(B, rep(0, m))
### makeH produces output A such that Ax >= 0
#print(qux)
out <- scdd(d2q(qux),adjacency = T,incidence = T,inputincidence = T,inputadjacency = T, representation = "H")
out1 = q2d(out$output)
#Why removing first two columns? Is it specific to this simulation setting? Or should it be for any A?

Delta = t(out1[,-c(1:2)])
dim(Delta)

adjacency = out$adjacency
adj.list = graph_from_adj_list(adjacency, mode = c("out", "in", "all", "total"),duplicate = TRUE)
adj.mat = as_adjacency_matrix(adj.list, sparse=F)

### cliques 
graph.list = graph_from_adjacency_matrix(adj.mat, "undirected", weighted=T, diag=F)
clique.list = max_cliques(graph.list)
clique.size = unlist(lapply(clique.list,length))
degree = clique.size

l = count_max_cliques(graph.list)

gamma_ls <- list()
Ti=24
#n = nrow(ydata3.scale)
n = nrow(ydata1)
xt = rep(1:Ti,n)
y = as.numeric(t(ydata1))
#y = array(ydata3.scale)
BS <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
Deltatilde <- BS %*% Delta

Als <- list()

for(i in 1:n){
  #i=1
  ni <- 1:Ti
  Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
}

X <- Deltatilde
cliques=clique.list
Total_itr=30

id = rep(1:n,each=Ti)

out4 <- fitting(y, X, BS1, cliques=clique.list, Total_itr=100, id, Als)
out2 <- fittingUN(y, xt, Ti, Total_itr=30)

out1$esti
out2$esti

par(mfrow=c(1,1))
plot(out1$gamma)
plot(out1$esti,type="l",ylim=c(0,355))
points(colMeans(ydata1))

plot(out2$esti,type="l",ylim=c(0,355))
points(colMeans(ydata1))
