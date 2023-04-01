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

source("/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava/fittingCompoSpikeFinal - Mixed effect.R")
source('/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava/fittingCompo - Mixed effect.R')

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

# b.ind = as.vector(clique.list[[2]])
# ftf = rowMeans(Delta[,b.ind])
plot(x,ftf,type="b",ylab="f(x)")
#plot(x,y,type="b",ylab="f(x)")
check = A%*%ftf


###############model fitting
J    <- 20
knot <- J-3

BS1 <- bsplineS((1:Ti)/Ti, breaks = seq(0, 1, 1/knot))

#A <- A %*% BS1

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
rep0 = 50

err <- array(0, c(rep0, 2, 3))
esti1 <- array(0, c(rep0, Ti, 3))
esti2 <- array(0, c(rep0, Ti, 3))
sigl <- c(0.5, 1, 2)
for(j in 1:3){
  print(j)
  for(rep in 1:rep0)
  {
    print(rep)
    #Generate data
    #j=1;rep=1
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
    
    #BS <- matrix(0, length(y), Ti)
    #BS[cbind(1:length(y), xt)] <- 1#bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
    BS <- bsplineS(xt/Ti, breaks = seq(0, 1, 1/knot))
    Deltatilde <- BS %*% Delta
    
    Als <- list()
    
    for(i in 1:n){
      ni <- which(id==i)
      
      Als[[i]] <- bsplineS(xt[ni]/Ti, breaks = seq(0, 1, 1/knot))
    }
    
    X <- Deltatilde
    cliques=clique.list
    Total_itr=30
    
    out2 <- fittingUN(y, xt, Ti, Total_itr=30)
    out1 <- fitting(y, X, BS1, cliques=clique.list, Total_itr=30)
    
    err[rep, 1, j] <- out1$err
    err[rep, 2, j] <- out2$err
    
    esti1[rep, , j] <- out1$esti
    esti2[rep, , j] <- out2$esti
    
    gamma_ls[[rep+(j-1)*rep0]] <- out1$gamma
    #print(err[rep, ])
  }
}

### 1st row spline + restricted
### 2nd row spline + unrestricted
err.median = apply(err, 2:3, median)
err.mean = apply(err, 2:3, mean)
err = round(rbind(err.median,err.mean),3)
colnames(err) = c("sd=0.5","sd=1","sd=2")
rownames(err) = c("median RS","median URS","mean RS","mean URS")
err
write.csv(err,"Bellshaped_spline_dense.csv")

### esti1 columns are each time point; rows are reps; j=3 matrices for 3 diff sd values
estim <- apply(esti1, 2:3, median)
estim2 <- apply(esti2, 2:3, median)

library(reshape)
library(ggplot2)
dat.true = data.frame(x=1:Ti,y = ftf)
dat.sd1 = data.frame(x=1:Ti,RS=estim[,1],URS=estim2[,1])
dat.sd1 = melt(dat.sd1,id.vars = "x")
dat.sd1$sd = "sd=0.5"
dat.sd2 = data.frame(x=1:Ti,RS=estim[,2],URS=estim2[,2])
dat.sd2 = melt(dat.sd2,id.vars = "x")
dat.sd2$sd = "sd=1"
dat.sd3 = data.frame(x=1:Ti,RS=estim[,3],URS=estim2[,3])
dat.sd3 = melt(dat.sd3,id.vars = "x")
dat.sd3$sd = "sd=2"
dat = rbind(dat.sd1,dat.sd2,dat.sd3)
p1 = ggplot()+
  geom_line(data = dat,aes(x=x,y=value,color=variable))+facet_wrap(~sd)+
  #  scale_linetype_manual(values=c("dashed","solid","solid","solid"))+
  theme(text = element_text(size = 35))+
  geom_point(data = dat.true,aes(x=x,y=y))
#theme_classic()+
#theme(legend.position=c(0.88,0.15),legend.key.size = unit(1.5, 'cm'),
#      legend.title = element_blank())
ggsave('Bellshaped_spline_dense.png', p1, width = 17,height=11, dpi = 300)


plot(gamma_ls[[1]])
plot(gamma_ls[[7]])

which.max(gamma_ls[[9]])
