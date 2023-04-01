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
y <- as.matrix(PAXINTEN_C[, 6:ncol(PAXINTEN_C)])
ydata <- apply(y, 1, timevec)
ydata <- t(ydata)

### 7176 participants for 7 days = 50232
subfreq <- array(table(PAXINTEN_C$SEQN))
sub  <- unique(PAXINTEN_C$SEQN)

covar_ls <- rnhanesdata::process_covar(waves=c("C","D"),
                                       varnames = c("SDDSRVYR","WTMEC2YR", "WTINT2YR",
                                                    "SDMVPSU", "SDMVSTRA",
                                                    "RIDAGEMN", "RIDAGEEX", "RIDAGEYR", "RIDRETH1", "RIAGENDR",
                                                    "BMXWT", "BMXHT", "BMXBMI", "DMDEDUC2",
                                                    "ALQ101", "ALQ110", "ALQ120Q","ALQ120U", "ALQ130", "SMQ020", "SMD030", "SMQ040",
                                                    "MCQ220","MCQ160F", "MCQ160B", "MCQ160C",
                                                    "PFQ049","PFQ054","PFQ057","PFQ059", "PFQ061B", "PFQ061C", "DIQ010"))

#y <- apply(PAXINTEN_C[, 6:ncol(PAXINTEN_C)], 1, mean)

Covariate_C <- covar_ls$Covariate_C
Covariate_D <- covar_ls$Covariate_D

indr <- which(Covariate_C$SEQN %in% PAXINTEN_C$SEQN)
Covariate_Cr <- Covariate_C[indr, ]

Covariate_Cr$EducationAdult <- factor(Covariate_Cr$DMDEDUC2, levels=c(1,2,3,4,5,7,9),
                                     labels=c("Less than 9th grade","9-11th grade","High school grad/GED or equivalent",
                                              "Some College or AA degree", "College graduate or above","Refused","Don't know"), ordered=FALSE)

#Inclusion criteria
ind1 <- which(as.numeric(Covariate_Cr$RIDAGEYR) <40 & as.numeric(Covariate_Cr$RIDAGEYR) > 30 & Covariate_Cr$EducationAdult == "College graduate or above")

samsub <- Covariate_Cr$SEQN[ind1]

#samsub <- sample(samsub, 200)

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

ydata2 <- ydatar[Sind2, ]
ind1 <- which(is.na(rowMeans(ydata2)==T))
if(length(ind1)) {ydata2 <- ydata2[-ind1, ]}

ydata3 <- ydatar[Sind3, ]
ind1 <- which(is.na(rowMeans(ydata3)==T))
if(length(ind1)) {ydata3 <- ydata1[-ind1, ]}

ydata4 <- ydatar[Sind4, ]
ind1 <- which(is.na(rowMeans(ydata4)==T))
if(length(ind1)) {ydata4 <- ydata4[-ind1, ]}

ydata5 <- ydatar[Sind5, ]
ind1 <- which(is.na(rowMeans(ydata5)==T))
if(length(ind1)) {ydata5 <- ydata5[-ind1, ]}

ydata6 <- ydatar[Sind6, ]
ind1 <- which(is.na(rowMeans(ydata6)==T))
if(length(ind1)) {ydata6 <- ydata6[-ind1, ]}

ydata7 <- ydatar[Sind7, ]
ind1 <- which(is.na(rowMeans(ydata7)==T))
if(length(ind1)) {ydata7 <- ydata7[-ind1, ]}

#Avg. activity over the week
ydataF <- (ydata1+ydata2+ydata3+ydata4+ydata5+ydata6+ydata7)/7
#range(as.numeric(Covariate_Cr$RIDAGEYR)[ind1])

#_C and _D correspond to the 2003-2004 and 2005-2006 
#
# SEQN: Respondent sequence number
# 
# WEEKDAY: Day of the week; WEEKDAY=1 for Sunday, 2 for Monday and so forth
# 
# PAXCAL: Denotes whether the monitor was in calibration when it was returned by the subject. The data for monitors that were out of calibration (PAXCAL=2) may be less reliable.
# 
# PAXSTAT: Component status code with PAXSTAT=1 for records with data that are deemed reliable. A PAXSTAT=2 was used to code records that had some questionable data; analysts may wish to examine these records more closely.
# 
# SDDSRVYR: Variable indicating which wave of the NHANES study this data is associated with. For example, SDDSRVYR = 3 corresponds to the 2003-2004 wave and SDDSRVYR = 4 corresponds to the 2005-2006 wave.


# ind <- which(PAXINTEN_C$SEQN == Covariate_Cr$SEQN[ind1[3]])
# plot(x[ind],y[ind])
# 
# 
# library(tidyverse)
# library(hrbrthemes)
# library(kableExtra)
# options(knitr.table.format = "html")
# #library(babynames)
# #library(streamgraph)
# library(viridis)
# library(DT)
# library(plotly)
# 
# v <- PAXINTEN_C[2000, 6:ncol(PAXINTEN_C)]
# 
# plot(log(1+unlist(v)), type='l')
# 
# v <- PAXINTEN_D[2000, 6:ncol(PAXINTEN_C)]
# 
# plot(log(1+unlist(v)), type='l')
# 
# y <- apply(PAXINTEN_D[, 6:ncol(PAXINTEN_D)], 1, function(x){diff(range(x))})
# x <- PAXINTEN_D[, 4]
# 
# ind <- which(PAXINTEN_D$SEQN == sub1[900])
# 
# plot(x[ind],y[ind])
# 
# 
# mat <- cbind(log(1+as.numeric(y)), as.numeric(x), as.numeric(PAXINTEN_D$SEQN))
# 
# colnames(mat) <- c("active", "day", "subj")
# data <- data.frame(mat)
# 
# data %>%
#   ggplot( aes(x=day, y=active, group=subj, color=subj)) +
#   geom_line() +
#   scale_color_viridis(discrete = F) +
#   theme(
#     legend.position="none",
#     plot.title = element_text(size=14)
#   ) +
#   ggtitle("A spaghetti chart of CN") #+
# #theme_ipsum()