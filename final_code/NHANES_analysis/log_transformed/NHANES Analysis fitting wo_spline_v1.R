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

Ti=24
t=1:24
i1 = 8
A = matrix(0,nrow=Ti+2,ncol=Ti)
for(j in 1:(i1-2)){
  A[j,j] = t[(j+2)]-t[(j+1)]
  A[j,(j+1)] = t[j]-t[(j+2)]
  A[j,(j+2)] = t[j+1]-t[j]
}

for(j in (i1-1):(Ti-2)){
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

J    <- 24
knot <- J-3

gamma_ls <- list()
Ti=24
n = nrow(ydata1.log)
xt = rep(1:Ti,n)
y = as.numeric(t(ydata1.log))
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
### ma and wo ma
out1 <- fitting(y, X, diag(Ti), cliques=clique.list, Total_itr=30, id, Als)
out2 <- fittingUN(y, xt, Ti, Total_itr=30)
out3 <- fittingUN(y, xt, Ti, splin=F, Total_itr=30)

result1 = data.frame(x=1:Ti, R = exp(out1$esti),
                    URS = exp(out2$esti), UR = exp(out3$esti))
result = cbind(result1,exp(out4$esti))
colnames(result)[5] = "RS"
#write.csv(result , "result_day1_exp.csv")

dat.true = data.frame(x=1:Ti,y=exp(colMeans(ydata1.log)))
df.plot = melt(result,id.vars = "x")
p1 = ggplot()+
  geom_line(data = df.plot,aes(x=x,y=value,color=variable))+
  #  scale_linetype_manual(values=c("dashed","solid","solid","solid"))+
  theme(text = element_text(size = 30))+
  geom_point(data = dat.true,aes(x=x,y=y))+ ylab("activity") + xlab("Time (in hours)")
p1  
ggsave('exp_log_ydata1.png', p1, width = 17,height=11, dpi = 300)

# par(mfrow=c(1,1))
# plot(out1$gamma)
# par(mfrow=c(2,1))
plot(out1$esti,type="l")
points(colMeans(ydata1.log))
# 
plot(out2$esti,type="l")
points(colMeans(ydata1.log))
# 
plot(out3$esti,type="l",ylim=c(0,5))
points(colMeans(ydata1.log))
# 
# out1$esti
# out2$esti
# out3$esti
# 
out1$sigma
out2$sigma
out3$sigma
