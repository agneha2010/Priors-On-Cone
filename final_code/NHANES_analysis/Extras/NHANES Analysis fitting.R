library(rnhanesdata)
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

temp = ydata[-ind,]
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
subfreq <- array(table(PAXINTEN_D$SEQN))
sub  <- unique(PAXINTEN_D$SEQN)

### covariates for defining inclusion criteria
indr <- which(Covariate_D$SEQN %in% PAXINTEN_D$SEQN)
Covariate_Dr <- Covariate_D[indr, ]

#Inclusion criteria
ind1 <- which(as.numeric(Covariate_Dr$RIDAGEYR) < 45 & 
                as.numeric(Covariate_Dr$RIDAGEYR) > 30 & 
                Covariate_Dr$EducationAdult == "College graduate or above" &
                #Covariate_Dr$DrinkStatus == "Non-Drinker"
                Covariate_Dr$DrinkStatus == "Moderate Drinker"
)
length(ind1)
samsub <- Covariate_Dr$SEQN[ind1]

### back to activity data
ind <- which(PAXINTEN_D$SEQN[-ind] %in% samsub)
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

source("/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava/fittingCompoSpikeFinal - Mixed effect Data.R")
source('/Users/pro/Desktop/Backup/priors_on_cone/projects/arkaprava/fittingCompo - Mixed effect Data.R')


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
n = nrow(ydata3.scale)
xt = rep(1:Ti,n)
y = array(ydata1.log)
y = array(ydata3.scale)
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

out1 <- fitting(y, X, diag(Ti), cliques=clique.list, Total_itr=30)
out2 <- fittingUN(y, xt, Ti, Total_itr=30)
out3 <- fittingUN(y, xt, Ti, splin=F, Total_itr=30)

out1$err
out2$err
out3$err

out1$esti
out2$esti
out3$esti

par(mfrow=c(1,1))
plot(out1$gamma)
plot(out1[[2]])
plot(out2)
plot(out3)

## need to modify the plot for data analysis
### err columns are each approach; rows are reps; j=3 matrices for 3 diff sd values
### 1st row restricted
### 2nd row spline + unrestricted
### 3rd row unrestricted + no spline 
err.median = apply(err, 2:3, median)
err.mean = apply(err, 2:3, mean)
err = round(rbind(err.median,err.mean),3)
colnames(err) = c("sd=0.5","sd=1","sd=2")
rownames(err) = c("median R","median URS","median UR","mean R","mean URS","mean UR")
err
write.csv(err,"Bellshaped_wo_spline_clique.csv")

### esti1 columns are each time point; rows are reps; j=3 matrices for 3 diff sd values
estim <- apply(esti1, 2:3, median)
estim2 <- apply(esti2, 2:3, median)
estim3 <- apply(esti3, 2:3, median)

library(reshape)
library(ggplot2)
dat.true = data.frame(x=1:Ti,y = ftf)
dat.sd1 = data.frame(x=1:Ti,R=estim[,1],URS=estim2[,1],UR=estim3[,1])
dat.sd1 = melt(dat.sd1,id.vars = "x")
dat.sd1$sd = "sd=0.5"
dat.sd2 = data.frame(x=1:Ti,R=estim[,2],URS=estim2[,2],UR=estim3[,2])
dat.sd2 = melt(dat.sd2,id.vars = "x")
dat.sd2$sd = "sd=1"
dat.sd3 = data.frame(x=1:Ti,R=estim[,3],URS=estim2[,3],UR=estim3[,3])
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
ggsave('Bellshaped_wo_spline_clique.png', p1, width = 17,height=11, dpi = 300)


plot(gamma_ls[[1]])
plot(gamma_ls[[7]])

which.max(gamma_ls[[9]])




