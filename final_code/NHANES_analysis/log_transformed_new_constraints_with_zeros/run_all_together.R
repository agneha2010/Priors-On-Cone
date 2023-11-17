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
library(dplyr)
library(reshape2)
library(ggplot2)

###############model fitting

source("/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava_v1/fittingCompoSpikeFinal - Mixed effect Data.R")
source('/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava_v1/fittingCompo - Mixed effect Data.R')

ydatar = read.csv("ydatar.csv")
ydatar = ydatar[,-1]

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
day1 = apply(ydata1,2,median)
day1.mean = colMeans(ydata1)
plot(day1.mean,type="p")
plot(day1,type="p")

ydata1.log = log(ydata1+1)

#### setting constraints
Ti=24
n = nrow(ydata1.log)
y = as.numeric(t(ydata1.log))
xt = rep(1:Ti,n)

t=1:24
i1 = 5
i2 = 8
i3 = 20
A = matrix(0,nrow=Ti+2,ncol=Ti)

for(j in 1:i1){
  A[j,j] = 1.0
}

A[(i1+1),i1] = -1.0
A[(i1+1),(i1+1)] = 1.0

for(j in i1:(i2-1)){
  A[j+2,j] = t[(j+2)]-t[(j+1)]
  A[j+2,(j+1)] = t[j]-t[(j+2)]
  A[j+2,(j+2)] = t[j+1]-t[j]
}

for(j in (i3):(Ti-2)){
  A[j+2,j] = t[(j+2)]-t[(j+1)]
  A[j+2,(j+1)] = t[j]-t[(j+2)]
  A[j+2,(j+2)] = t[j+1]-t[j]
}

for(j in (i2):(i3-1)){
  A[j+2,j] = t[(j+1)]-t[(j+2)]
  A[j+2,(j+1)] = -t[j]+t[(j+2)]
  A[j+2,(j+2)] = -t[j+1]+t[j]
}

A[(Ti+1),(Ti-1)] = 1.0
A[(Ti+1),Ti] = -1.0
A[(Ti+2),Ti] = 1.0

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

J    <- 24
knot <- J-3

BS <- matrix(0, length(y), Ti)
BS[cbind(1:length(y), xt)] <- 1
Deltatilde <- BS %*% Delta

gamma_ls <- list()
Als <- list()

for(i in 1:n){
  ni <- 1:Ti
  Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
}

X <- Deltatilde
cliques=clique.list
Total_itr=30

id = rep(1:n,each=Ti)
### ma and wo ma
R <- fitting(y, X, diag(Ti), cliques=clique.list, Total_itr=100, id, Als)
URS <- fittingUN(y, xt, Ti, Total_itr=100)
UR <- fittingUN(y, xt, Ti, splin=F, Total_itr=100)

result = data.frame(x=1:Ti,R = R$esti,URS = URS$esti,UR = UR$esti)

#### spline part

t = 1:24
i1 = 5
i2 = 8
i3 = 20
A = matrix(0,nrow=Ti+2,ncol=Ti)

for(j in 1:i1){
  A[j,j] = 1.0
}

A[(i1+1),i1] = -1.0
A[(i1+1),(i1+1)] = 1.0

for(j in i1:(i2-1)){
  A[j+2,j] = t[(j+2)]-t[(j+1)]
  A[j+2,(j+1)] = t[j]-t[(j+2)]
  A[j+2,(j+2)] = t[j+1]-t[j]
}

for(j in (i3):(Ti-2)){
  A[j+2,j] = t[(j+2)]-t[(j+1)]
  A[j+2,(j+1)] = t[j]-t[(j+2)]
  A[j+2,(j+2)] = t[j+1]-t[j]
}

for(j in (i2):(i3-1)){
  A[j+2,j] = t[(j+1)]-t[(j+2)]
  A[j+2,(j+1)] = -t[j]+t[(j+2)]
  A[j+2,(j+2)] = -t[j+1]+t[j]
}

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

BS <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
Deltatilde <- BS %*% Delta

gamma_ls <- list()
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

RS <- fitting(y, X, BS1, cliques=clique.list, Total_itr=100, id, Als)
#URS1 <- fittingUN(y, xt, Ti, Total_itr=30)

result = cbind(result,RS$esti)
colnames(result)[5] = c("RS")
write.csv(result , "ydata7.log.csv")


dat.true = data.frame(x=1:Ti,y=colMeans(ydata7.log))
df.plot = melt(result,id.vars = "x")
p1 = ggplot()+
  geom_line(data = df.plot,aes(x=x,y=value,color=variable))+
  #  scale_linetype_manual(values=c("dashed","solid","solid","solid"))+
  theme(text = element_text(size = 30))+
  geom_point(data = dat.true,aes(x=x,y=y))+ ylab("activity") + xlab("Time (in hours)")
p1  
ggsave('ydata7.log.png', p1, width = 17,height=11, dpi = 300)

