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

## load the data
data("PAXINTEN_C");data("PAXINTEN_D")
data("Flags_C");data("Flags_D")
#data("Mortality_2011_C");data("Mortality_2011_D")
data("Covariate_C");data("Covariate_D")
data("Mortality_2015_C");data("Mortality_2015_D")

## 60*24=1440 minute level data to 24 hour data
timevec <- function(x){
  mat <- matrix(x, 60, 24)
  
  return(array(colMeans(mat, na.rm = T)))
}
Flags_D = as.matrix(Flags_D)
Flags_Dr = ifelse(Flags_D[, 6:ncol(Flags_D)]==0,NA,Flags_D[, 6:ncol(Flags_D)])
y <- as.matrix(PAXINTEN_D[, 6:ncol(PAXINTEN_D)]*Flags_Dr)
ydata <- apply(y, 1, timevec)
ydata <- t(ydata)
ind = which(is.na(rowMeans(ydata)))

ydata = ydata[-ind,]
for(i in 1:6){
  plot(temp[i,])
}

#timevec.med <- function(x){
#   mat <- matrix(x, 60, 24)
#   
#   return(array(apply(mat,2,median,na.rm=T)))
# }
#ydata.med <- apply(y, 1, timevec.med)
# ydata.med <- t(ydata.med)
# ind.med = which(is.na(rowMeans(ydata.med)))
# 
# temp.med = ydata.med[-ind,]
# for(i in 1:3){
#   plot(temp[i,])
#   plot(temp.med[i,])
# }

### 7176 participants for 7 days = 50232
ydatar <- ydata

Sind1 <- (0:(nrow(ydatar)/7-1))*7+1
Sind2 <- (0:(nrow(ydatar)/7-1))*7+2
Sind3 <- (0:(nrow(ydatar)/7-1))*7+3
Sind4 <- (0:(nrow(ydatar)/7-1))*7+4
Sind5 <- (0:(nrow(ydatar)/7-1))*7+5
Sind6 <- (0:(nrow(ydatar)/7-1))*7+6
Sind7 <- (0:(nrow(ydatar)/7-1))*7+7

ydata1 <- ydatar[Sind1, ]
ind1 <- which(is.na(rowMeans(ydata1)==T))
if(length(ind1)) {ydata1 <- ydata1[-ind1, ]}

ydata3 <- ydatar[Sind3, ]
ind3 <- which(is.na(rowMeans(ydata3)==T))
if(length(ind3)) {ydata3 <- ydata3[-ind3, ]}

ydata1.log = log(ydata1+1)
ydata1.scale = t(apply(ydata1,1,function(x) (x-min(x))/(max(x)-min(x))))
apply(ydata1.log,2,sd)

ydata3.log = log(ydata3+1)
ydata3.scale = t(apply(ydata3,1,function(x) (x-min(x))/(max(x)-min(x))))
apply(ydata3.log,2,sd)
ind = which(is.na(rowMeans(ydata3.scale)))
ydata3.scale = ydata3.scale[-ind,]
###############model fitting

source("/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava_v1/fittingCompoSpikeFinal - Mixed effect Data.R")
source('/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava_v1/fittingCompo - Mixed effect Data.R')

m = nrow(A)
### makeH uses Ax <= 0; so multiply by -1
B = -A
qux <- makeH(B, rep(0, m))
### makeH produces output A such that Ax >= 0
#print(qux)
out <- scdd(d2q(qux),adjacency = T,incidence = T,inputincidence = T,inputadjacency = T, representation = "H")
out1 = q2d(out$output)

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

# b.ind = as.vector(clique.list[[2]])
# ftf = rowMeans(Delta[,b.ind])
# plot(x,ftf,type="b",ylab="f(x)")

J    <- 24
knot <- J-3

gamma_ls <- list()
Ti=24
#n = nrow(ydata3.scale)
n = nrow(ydata3.log)
xt = rep(1:Ti,n)
y = array(ydata3.log)
#y = array(ydata3.scale)
BS <- matrix(0, length(y), Ti)
BS[cbind(1:length(y), xt)] <- 1
Deltatilde <- BS %*% Delta

Als <- list()

for(i in 1:n){
  ni <- 1:Ti
  
  Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
}

X <- Deltatilde
cliques=clique.list
Total_itr=30

id = rep(1:n,each=Ti)

out1 <- fitting(y, X, diag(Ti), cliques=clique.list, Total_itr=30, id, Als)
out2 <- fittingUN(y, xt, Ti, Total_itr=30)
out3 <- fittingUN(y, xt, Ti, splin=F, Total_itr=30)


out1$esti
out2$esti
out3$esti

par(mfrow=c(1,1))
plot(out1$gamma)
plot(out1$esti)
plot(out2$esti)
plot(out3$esti)
plot(ydata3.log[1,])
