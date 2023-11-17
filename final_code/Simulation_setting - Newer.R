# #PAXINTEN_D <- read.csv("~/NHANES/PAXINTEN_D.csv", header = T)
# 
# # ind <- which(PAXINTEN_D$SEQN==31193)
# # 
# # y <- as.matrix(PAXINTEN_D[ind, 6:ncol(PAXINTEN_D)])
# # x <- 1:1440#PAXINTEN_C[, 4]
# 
# 
# require(stats)
# require(graphics)
# library(splines2)
# library(coneproj)
# library(rcdd)
# library(truncnorm)
# library(truncdist)
# library(mvtnorm)
# library(R.utils)
# library(Matrix)
# library(igraph)
# library(MCMCpack)
# library(fda)
# 
# source("~/Conic/fittingCompoSpikeFinal - Mixed effect.R")
# source("~/Conic/fittingCompo - Mixed effect.R")
# 
# n <- 100
# Ti <- 20
# Tiv <- 1:20
# 
# sd = 0.5;m=0
# x = seq(-2,2,length.out = Ti)
# x = sort(x)
# set.seed(11131991)
# fx = dnorm(x,mean=m,sd=sd) ### true mean
# plot(x,fx)
# 
# ftf <- 10*fx
# 
# ## finding the A matrix for bell-shaped
# d2fx = function(x){
#   return(dnorm(x,mean=m,sd=sd)*(-1+(x-m)^2/sd^2)/sd^2)
# }
# d2fx(x)
# ind = which(d2fx(x)<0)
# t=x
# i1 = ind[1]
# i2 = ind[length(ind)]
# # 20 30
# A = matrix(0,nrow=Ti+2,ncol=Ti)
# for(j in 1:(i1-2)){
#   A[j,j] = t[(j+2)]-t[(j+1)]
#   A[j,(j+1)] = t[j]-t[(j+2)]
#   A[j,(j+2)] = t[j+1]-t[j]
# }
# 
# for(j in (i2):(Ti-2)){
#   A[j,j] = t[(j+2)]-t[(j+1)]
#   A[j,(j+1)] = t[j]-t[(j+2)]
#   A[j,(j+2)] = t[j+1]-t[j]
# }
# 
# for(j in (i1-1):(i2-1)){
#   A[j,j] = t[(j+1)]-t[(j+2)]
#   A[j,(j+1)] = -t[j]+t[(j+2)]
#   A[j,(j+2)] = -t[j+1]+t[j]
# }
# ## added by me
# #A[i2-1,(i2-2):i2] = - A[i2-1,(i2-2):i2]
# 
# A[(Ti-1),1] = -1.0
# A[(Ti-1),2] = 1.0
# A[Ti,1] = 1.0
# A[(Ti+1),(Ti-1)] = 1.0
# A[(Ti+1),Ti] = -1.0
# A[(Ti+2),Ti] = 1.0
# 
# 
# ###############model fitting
# J    <- 20
# knot <- J-3
# 
# BS1 <- bsplineS((1:Ti)/Ti, breaks = seq(0, 1, 1/knot))
# 
# A <- A %*% BS1
# 
# m = nrow(A)
# ### makeH uses Ax <= 0; so multiply by -1
# B = -A
# qux <- makeH(B, rep(0, m))
# ### makeH produces output A such that Ax >= 0
# #print(qux)
# out <- scdd(d2q(qux),adjacency = T,incidence = T,inputincidence = T,inputadjacency = T, representation = "H")
# out1 = q2d(out$output)
# #Why removing first two columns? Is it specific to this simulation setting? Or should it be for any A?
# 
# Delta = t(out1[,-c(1:2)])
# dim(Delta)
# 
# 
# adjacency = out$adjacency
# adj.list = graph_from_adj_list(adjacency, mode = c("out", "in", "all", "total"),duplicate = TRUE)
# adj.mat = as_adjacency_matrix(adj.list, sparse=F)
# 
# ### cliques 
# graph.list = graph_from_adjacency_matrix(adj.mat, "undirected", weighted=T, diag=F)
# clique.list = max_cliques(graph.list)
# clique.size = unlist(lapply(clique.list,length))
# degree = clique.size
# 
# l = count_max_cliques(graph.list)
# 
# #ft[7:Ti] <- ft[7]
# 
# #plot(ft)
# 
# gamma_ls <- list()
# 
# err <- array(0, c(50, 2, 3))
# esti1 <- array(0, c(50, Ti, 3))
# esti2 <- array(0, c(50, Ti, 3))
# sigl <- c(1, 5, 10)
# for(j in 1:3){
#   for(rep in 1:50)
#   {
#     #Generate data
#     
#     xt <- rep(1:Ti, n)
#     id <- 1
#     # while(length(table(id)) < n){
#     #   ind <- sample(1:(Ti*n), Ti*n)
#     #   id <- rep(1:n, each=Ti)
#     #   id <- id[-ind]
#     # }
#     id <- rep(1:n, each=Ti)
#     xt <- xt#[-ind]
#     
#     fxt <- ftf[xt]#exp(-abs(xt-2.6)/3) + a^(-abs(xt-2.6)/3)
#     
#     #fxt[which(xt <=10 & xt>=7)] <- ft[7]
#     
#     #plot(fxt)
#     
#     Sigma <- exp(-as.matrix(dist(1:Ti))/4)
#     
#     W <- mvtnorm::rmvnorm(n, sigma = Sigma)
#     W <- array(t(W))#[-ind]
#     
#     sige <- sigl[j]
#     e <- rnorm(length(xt), sd=sige)
#     
#     y <- array(fxt + W + e)
#     
#     #A <- rbind(A, diag(J))
#     
#     BS <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
#     Deltatilde <- BS %*% Delta
#     
#     Als <- list()
#     
#     for(i in 1:n){
#       ni <- which(id==i)
#       
#       Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
#     }
#     
#     X <- Deltatilde
#     cliques=clique.list
#     Total_itr=30
#     
#     out2 <- fittingUN(y, xt, Ti, Total_itr=30)
#     out1 <- fitting(y, X, BS1, cliques=clique.list, Total_itr=30)
#     
#     err[rep, 1, j] <- out1$err
#     err[rep, 2, j] <- out2$err
#     
#     esti1[rep, , j] <- out1$esti
#     esti2[rep, , j] <- out2$esti
#     
#     gamma_ls[[rep+(j-1)*50]] <- out1$gamma
#     #print(err[rep, ])
#   }
# }
# 
# apply(err, 2:3, median)
# apply(err, 2:3, mean)
# setwd("~/Conic")
# save(err, esti1, esti2, gamma_ls, file="result10F.rda")
# 
# estim <- apply(esti1, 2:3, mean)
# plot(1:Ti, ftf, type='l', ylab="f(x)", xlab="Time", cex.lab=1.5)
# points(estim[, 1], type='l', col=2)
# points(estim[, 2], type='l', col=3)
# points(estim[, 3], type='l', col=4)
# legend("topleft",c("True f(x)", "Error SD = 1", "Error SD = 5", "Error SD = 10"),lwd=2, ncol=1, col=1:4)
# #######################################################################################
# n <- 100
# Ti <- 20
# Tiv <- 1:20
# 
# sd = 0.5;m=0
# x = seq(-2,2,length.out = Ti)
# x = sort(x)
# set.seed(11131991)
# fx = dnorm(x,mean=m,sd=sd) ### true mean
# plot(x,fx)
# 
# ftf <- 3*fx
# 
# ## finding the A matrix for bell-shaped
# d2fx = function(x){
#   return(dnorm(x,mean=m,sd=sd)*(-1+(x-m)^2/sd^2)/sd^2)
# }
# d2fx(x)
# ind = which(d2fx(x)<0)
# t=x
# i1 = ind[1]
# i2 = ind[length(ind)]
# # 20 30
# A = matrix(0,nrow=Ti+2,ncol=Ti)
# for(j in 1:(i1-2)){
#   A[j,j] = t[(j+2)]-t[(j+1)]
#   A[j,(j+1)] = t[j]-t[(j+2)]
#   A[j,(j+2)] = t[j+1]-t[j]
# }
# 
# for(j in (i2):(Ti-2)){
#   A[j,j] = t[(j+2)]-t[(j+1)]
#   A[j,(j+1)] = t[j]-t[(j+2)]
#   A[j,(j+2)] = t[j+1]-t[j]
# }
# 
# for(j in (i1-1):(i2-1)){
#   A[j,j] = t[(j+1)]-t[(j+2)]
#   A[j,(j+1)] = -t[j]+t[(j+2)]
#   A[j,(j+2)] = -t[j+1]+t[j]
# }
# ## added by me
# #A[i2-1,(i2-2):i2] = - A[i2-1,(i2-2):i2]
# 
# A[(Ti-1),1] = -1.0
# A[(Ti-1),2] = 1.0
# A[Ti,1] = 1.0
# A[(Ti+1),(Ti-1)] = 1.0
# A[(Ti+1),Ti] = -1.0
# A[(Ti+2),Ti] = 1.0
# 
# 
# ###############model fitting
# J    <- 20
# knot <- J-3
# 
# BS1 <- bsplineS((1:Ti)/Ti, breaks = seq(0, 1, 1/knot))
# 
# #A <- A %*% BS1
# 
# m = nrow(A)
# ### makeH uses Ax <= 0; so multiply by -1
# B = -A
# qux <- makeH(B, rep(0, m))
# ### makeH produces output A such that Ax >= 0
# #print(qux)
# out <- scdd(d2q(qux),adjacency = T,incidence = T,inputincidence = T,inputadjacency = T, representation = "H")
# out1 = q2d(out$output)
# #Why removing first two columns? Is it specific to this simulation setting? Or should it be for any A?
# 
# Delta = t(out1[,-c(1:2)])
# dim(Delta)
# 
# 
# adjacency = out$adjacency
# adj.list = graph_from_adj_list(adjacency, mode = c("out", "in", "all", "total"),duplicate = TRUE)
# adj.mat = as_adjacency_matrix(adj.list, sparse=F)
# 
# ### cliques 
# graph.list = graph_from_adjacency_matrix(adj.mat, "undirected", weighted=T, diag=F)
# clique.list = max_cliques(graph.list)
# clique.size = unlist(lapply(clique.list,length))
# degree = clique.size
# 
# l = count_max_cliques(graph.list)
# 
# b.ind = as.vector(clique.list[[2]])
# ftf = rowMeans(Delta[,b.ind])
# #ft[7:Ti] <- ft[7]
# 
# #plot(ft)
# 
# A <- A %*% BS1
# 
# m = nrow(A)
# ### makeH uses Ax <= 0; so multiply by -1
# B = -A
# qux <- makeH(B, rep(0, m))
# ### makeH produces output A such that Ax >= 0
# #print(qux)
# out <- scdd(d2q(qux),adjacency = T,incidence = T,inputincidence = T,inputadjacency = T, representation = "H")
# out1 = q2d(out$output)
# #Why removing first two columns? Is it specific to this simulation setting? Or should it be for any A?
# 
# Delta = t(out1[,-c(1:2)])
# dim(Delta)
# 
# 
# adjacency = out$adjacency
# adj.list = graph_from_adj_list(adjacency, mode = c("out", "in", "all", "total"),duplicate = TRUE)
# adj.mat = as_adjacency_matrix(adj.list, sparse=F)
# 
# ### cliques 
# graph.list = graph_from_adjacency_matrix(adj.mat, "undirected", weighted=T, diag=F)
# clique.list = max_cliques(graph.list)
# clique.size = unlist(lapply(clique.list,length))
# degree = clique.size
# 
# l = count_max_cliques(graph.list)
# 
# 
# gamma_ls <- list()
# 
# err <- array(0, c(50, 2, 3))
# esti1 <- array(0, c(50, Ti, 3))
# esti2 <- array(0, c(50, Ti, 3))
# sigl <- c(1, 2, 5)
# for(j in 1:3){
#   for(rep in 1:50)
#   {
#     #Generate data
#     
#     xt <- rep(1:Ti, n)
#     id <- 1
#     while(length(table(id)) < n){
#       ind <- sample(1:(Ti*n), 0.2*Ti*n)
#       id <- rep(1:n, each=Ti)
#       id <- id[-ind]
#     }
#     
#     xt <- xt[-ind]
#     
#     fxt <- ftf[xt]#exp(-abs(xt-2.6)/3) + a^(-abs(xt-2.6)/3)
#     
#     #fxt[which(xt <=10 & xt>=7)] <- ft[7]
#     
#     #plot(fxt)
#     
#     Sigma <- exp(-as.matrix(dist(1:Ti))/4)
#     
#     W <- mvtnorm::rmvnorm(n, sigma = Sigma)
#     W <- array(t(W))[-ind]
#     
#     sige <- sigl[j]
#     e <- rnorm(length(xt), sd=sige)
#     
#     y <- array(fxt + W + e)
#     
#     #A <- rbind(A, diag(J))
#     
#     BS <- matrix(0, length(y), Ti)
#     BS[cbind(1:length(y), xt)] <- 1#bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
#     BS <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
#     Deltatilde <- BS %*% Delta
#     
#     Als <- list()
#     
#     for(i in 1:n){
#       ni <- which(id==i)
#       
#       Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
#     }
#     
#     X <- Deltatilde
#     cliques=clique.list
#     Total_itr=30
#     
#     out2 <- fittingUN(y, xt, Ti, Total_itr=30)
#     out1 <- fitting(y, X, BS1, cliques=clique.list, Total_itr=30)
#     
#     err[rep, 1, j] <- out1$err
#     err[rep, 2, j] <- out2$err
#     
#     esti1[rep, , j] <- out1$esti
#     esti2[rep, , j] <- out2$esti
#     
#     gamma_ls[[rep+(j-1)*50]] <- out1$gamma
#     #print(err[rep, ])
#   }
# }
# 
# apply(err, 2:3, median)
# apply(err, 2:3, mean)
# setwd("~/Conic")
# save(err, esti1, esti2, gamma_ls, file="otherresultF.rda")
# 
# estim <- apply(esti1, 2:3, mean)
# plot(1:Ti, ftf, type='l', ylab="f(x)", xlab="Time", cex.lab=1.5)
# points(estim[, 1], type='l', col=2)
# points(estim[, 2], type='l', col=3)
# points(estim[, 3], type='l', col=4)
# legend("topleft",c("True f(x)", "Error SD = 1", "Error SD = 5", "Error SD = 10"),lwd=2, ncol=1, col=1:4)


#######################################################################################

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

source("~/Conic/fittingCompoSpikeFinal - Mixed effect.R")
source("~/Conic/fittingCompo - Mixed effect.R")


n <- 100
Ti <- 20
Tiv <- 1:20

sd = 0.5;m=0
x = seq(-2,2,length.out = Ti)
x = sort(x)
set.seed(11131991)
fx = dnorm(x,mean=m,sd=sd) ### true mean
plot(x,fx)

ftf <- 3*fx

## finding the A matrix for bell-shaped
d2fx = function(x){
  return(dnorm(x,mean=m,sd=sd)*(-1+(x-m)^2/sd^2)/sd^2)
}
d2fx(x)
ind = which(d2fx(x)<0)
t=x
i1 = ind[1]
i2 = ind[length(ind)]
# 20 30
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


###############model fitting
J    <- 20
knot <- J-3

BS1 <- bsplineS((1:Ti)/Ti, breaks = seq(0, 1, 1/knot))

#A <- A %*% BS1

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

b.ind = as.vector(clique.list[[2]])
ftf = rowMeans(Delta[,b.ind])
#ft[7:Ti] <- ft[7]

#plot(ft)

# A <- A %*% BS1
# 
# m = nrow(A)
# ### makeH uses Ax <= 0; so multiply by -1
# B = -A
# qux <- makeH(B, rep(0, m))
# ### makeH produces output A such that Ax >= 0
# #print(qux)
# out <- scdd(d2q(qux),adjacency = T,incidence = T,inputincidence = T,inputadjacency = T, representation = "H")
# out1 = q2d(out$output)
# #Why removing first two columns? Is it specific to this simulation setting? Or should it be for any A?
# 
# Delta = t(out1[,-c(1:2)])
# dim(Delta)
# 
# 
# adjacency = out$adjacency
# adj.list = graph_from_adj_list(adjacency, mode = c("out", "in", "all", "total"),duplicate = TRUE)
# adj.mat = as_adjacency_matrix(adj.list, sparse=F)
# 
# ### cliques 
# graph.list = graph_from_adjacency_matrix(adj.mat, "undirected", weighted=T, diag=F)
# clique.list = max_cliques(graph.list)
# clique.size = unlist(lapply(clique.list,length))
# degree = clique.size
# 
# l = count_max_cliques(graph.list)


gamma_ls <- list()

err <- array(0, c(50, 4, 3))
esti1 <- array(0, c(50, Ti, 3))
esti2 <- array(0, c(50, Ti, 3))
esti3 <- array(0, c(50, Ti, 3))
esti4 <- array(0, c(50, Ti, 3))
sigl <- c(1, 5, 10)
for(j in 1:1){
  for(rep in 1:50)
  {
    #Generate data
    
    xt <- rep(1:Ti, n)
    id <- 1
    while(length(table(id)) < n){
      ind <- sample(1:(Ti*n), 0.2*Ti*n)
      id <- rep(1:n, each=Ti)
      id <- id[-ind]
    }
    
    xt <- xt[-ind]
    
    fxt <- ftf[xt]#exp(-abs(xt-2.6)/3) + a^(-abs(xt-2.6)/3)
    
    #fxt[which(xt <=10 & xt>=7)] <- ft[7]
    
    #plot(fxt)
    
    Sigma <- exp(-as.matrix(dist(1:Ti))/4)
    
    W <- mvtnorm::rmvnorm(n, sigma = Sigma)
    W <- array(t(W))[-ind]
    
    sige <- sigl[j]
    e <- rnorm(length(xt), sd=sige)
    
    y <- array(fxt + W + e)
    
    #A <- rbind(A, diag(J))
    
    BS <- matrix(0, length(y), Ti)
    BS[cbind(1:length(y), xt)] <- 1#bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
    #BS <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
    Deltatilde <- BS %*% Delta
    
    Als <- list()
    
    for(i in 1:n){
      ni <- which(id==i)
      
      Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
    }
    
    X <- Deltatilde
    cliques=clique.list
    Total_itr=30
    
    out3 <- fittingUN(y, xt, Ti, Total_itr=Total_itr)
    out1 <- fitting(y, X, diag(Ti), cliques=clique.list, Total_itr=Total_itr, id, Als)
    #out2 <- fitting(y, X, BS1, cliques=clique.list, Total_itr=Total_itr)
    
    out4 <- fittingUN(y, xt, Ti, splin=F, Total_itr=Total_itr)
    
    err[rep, 1, j] <- out1$err
    err[rep, 2, j] <- out2$err
    err[rep, 3, j] <- out3$err
    err[rep, 4, j] <- out4$err
    
    esti1[rep, , j] <- out1$esti
    esti2[rep, , j] <- out2$esti
    esti3[rep, , j] <- out3$esti
    esti4[rep, , j] <- out4$esti
    
    gamma_ls[[rep+(j-1)*50]] <- out1$gamma
    #print(err[rep, ])
  }
}

apply(err, 2:3, median)
apply(err, 2:3, mean)
setwd("~/Conic")
save(err, esti1, esti2, gamma_ls, file="other1resultF.rda")

dev.new()
estim <- apply(esti1[1:15,,], 2:3, mean)
plot(1:Ti, ftf, type='l', ylab="f(x)", xlab="Time", cex.lab=1.5)
points(estim[, 1], type='l', col=2)
points(estim[, 2], type='l', col=3)
points(estim[, 3], type='l', col=4)
legend("topleft",c("True f(x)", "Error SD = 5", "Error SD = 10", "Error SD = 15"),lwd=2, ncol=1, col=1:4)

###################################################################
load("result10F.rda")

apply(err, 2:3, median)

plot(gamma_ls[[1]])
plot(gamma_ls[[21]])

load("otherresultF.rda")

apply(err, 2:3, median)

plot(gamma_ls[[1]])
plot(gamma_ls[[31]])


load("other1resultF.rda")

apply(err, 2:3, median)

plot(gamma_ls[[1]])
plot(gamma_ls[[71]])
