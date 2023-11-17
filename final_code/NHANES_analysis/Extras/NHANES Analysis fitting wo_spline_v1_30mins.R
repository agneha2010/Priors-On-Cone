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
library(zoo)

data("Flags_C");data("Flags_D")
#data("Mortality_2011_C");data("Mortality_2011_D")
data("Covariate_C");data("Covariate_D")
data("Mortality_2015_C");data("Mortality_2015_D")

## 60*24=1440 minute level data to 24 hour data
timevec <- function(x){
  mat <- matrix(x, 30, 48)
  
  return(array(colMeans(mat, na.rm = T)))
}
y <- as.matrix(PAXINTEN_C[, 6:ncol(PAXINTEN_C)])
ydata <- apply(y, 1, timevec)
ydata <- t(ydata)

### 7176 participants for 7 days = 50232
subfreq <- array(table(PAXINTEN_C$SEQN))
sub  <- unique(PAXINTEN_C$SEQN)

### covariates for defining inclusion criteria
indr <- which(Covariate_C$SEQN %in% PAXINTEN_C$SEQN)
Covariate_Cr <- Covariate_C[indr, ]

#Inclusion criteria
ind1 <- which(as.numeric(Covariate_Cr$RIDAGEYR) < 40 & 
                as.numeric(Covariate_Cr$RIDAGEYR) > 30 & 
                Covariate_Cr$EducationAdult == "College graduate or above" 
)
length(ind1)
samsub <- Covariate_Cr$SEQN[ind1]

### back to activity data
ind <- which(PAXINTEN_C$SEQN %in% samsub)
ydatar <- ydata[ind, ]

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

ydata3 <- ydatar[Sind2, ]
ind1 <- which(is.na(rowMeans(ydata3)==T))
if(length(ind1)) {ydata3 <- ydata3[-ind1, ]}

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

J    <- 48
knot <- J-3

gamma_ls <- list()
Ti=48
#n = nrow(ydata3.scale)
n = nrow(ydata2)
xt = rep(1:Ti,n)
y = array(ydata2)
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

plot(1:48,out1$esti,type="l",ylim=c(0,350))
points(1:48,colMeans(ydata2))

plot(out2$esti)
plot(out3$esti)
plot(ydata3.log[1,])
