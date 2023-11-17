library(rnhanesdata)
library(ggplot2)
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
day1 = apply(ydata1,2,median)
day1.mean = colMeans(ydata1)
plot(day1.mean,type="p")
plot(day1,type="p")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
ydata1 = as.data.frame(t(ydata1))
ydata.plot = ydata1[,1:20]
ydata.plot$x = 1:24
dat = melt(ydata.plot,id.vars = "x")
dat.mean = data.frame(x=1:24,y=day1.mean)
p1 = ggplot()+
  geom_line(data = dat,aes(x=x,y=value,color=variable),lwd=0.8)+
  geom_line(data = dat.mean,aes(x=x,y=y),size=1.5,linetype="dashed")+
  theme_classic()+
  theme(legend.position="None",legend.key.size = unit(1.5, 'cm'),
        legend.title = element_blank())+
  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))+
  theme(text = element_text(size = 30))+
  labs(x ="Time (in hours)", y = "Hourly activity Count")+
  scale_color_manual(values = c(brewer.pal(10, name="Spectral"),brewer.pal(10,  "Paired")))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.9))
p1
ggsave('20_subjects.png', p1, width = 17,height=11, dpi = 300)


ydata2 <- ydatar[Sind2, ]
ind1 <- which(is.na(rowMeans(ydata2)==T))
if(length(ind1)) {ydata2 <- ydata2[-ind1, ]}
day2 = apply(ydata2,2,median)
#day2 = colMeans(ydata2)
plot(day2,type="l")

ydata3 <- ydatar[Sind3, ]
ind1 <- which(is.na(rowMeans(ydata3)==T))
if(length(ind1)) {ydata3 <- ydata1[-ind1, ]}
day3 = apply(ydata3,2,median)
#day3 = colMeans(ydata3)
plot(day3,type="l")

ydata4 <- ydatar[Sind4, ]
ind1 <- which(is.na(rowMeans(ydata4)==T))
if(length(ind1)) {ydata4 <- ydata4[-ind1, ]}
day4 = apply(ydata4,2,median)
#day4 = colMeans(ydata4)
plot(day4,type="l")

ydata5 <- ydatar[Sind5, ]
ind1 <- which(is.na(rowMeans(ydata5)==T))
if(length(ind1)) {ydata5 <- ydata5[-ind1, ]}
day5 = apply(ydata5,2,median)
#day5 = colMeans(ydata5)
plot(day5,type="l")

ydata6 <- ydatar[Sind6, ]
ind1 <- which(is.na(rowMeans(ydata6)==T))
if(length(ind1)) {ydata6 <- ydata6[-ind1, ]}
day6 = apply(ydata6,2,median)
#day6 = colMeans(ydata6)
plot(day6,type="l")

ydata7 <- ydatar[Sind7, ]
ind1 <- which(is.na(rowMeans(ydata7)==T))
if(length(ind1)) {ydata7 <- ydata7[-ind1, ]}
day7 = apply(ydata7,2,median)
#day7 = colMeans(ydata7)
plot(day7,type="l")

par(mfrow=c(2,4))
plot(day1,type="l")
plot(day2,type="l")
plot(day3,type="l")
plot(day4,type="l")
plot(day5,type="l")
plot(day6,type="l")
plot(day7,type="l")

day1.med = data.frame(x=1:24,y=day1)
ggplot(data=day1.med,aes(x=x,y=y))+geom_point()+
  geom_smooth(method=loess, color="blue") 

day2.med = data.frame(x=1:24,y=day2)
ggplot(data=day2.med,aes(x=x,y=y))+geom_point()+
  geom_smooth(method=loess, color="blue")

day3.med = data.frame(x=1:24,y=day3)
ggplot(data=day3.med,aes(x=x,y=y))+geom_point()+
  geom_smooth(method=loess, color="blue")

day11 = data.frame(x=1:24,y=day1)
ggplot(data=day11,aes(x=x,y=y))+geom_point()+
  geom_smooth(method=loess, color="blue")


#check = A %*% day1.med[,2]
#View(cbind(check,A))
