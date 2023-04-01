## load the data from rnhanesdata package
data("PAXINTEN_C");data("PAXINTEN_D")
data("Flags_C");data("Flags_D")
data("Covariate_C");data("Covariate_D")
data("Mortality_2015_C");data("Mortality_2015_D")

## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact much data, most estimated non-wear times correspond to 0 counts anyway
PAXINTEN_C[,paste0("MIN",1:1440)] <- PAXINTEN_C[,paste0("MIN",1:1440)]*
  Flags_C[,paste0("MIN",1:1440)]
PAXINTEN_D[,paste0("MIN",1:1440)] <- PAXINTEN_D[,paste0("MIN",1:1440)]*
  Flags_D[,paste0("MIN",1:1440)]

## Merge covariate, mortality, and accelerometry data
## note that both PAXINTEN_* and Covariate_* have a column
## called "SDDSRVYR" indicating which NHANES wave the data is associated with.
## To avoid duplicating this column in the merged data, we add this variable to the "by"
## argument in left_join()
AllAct_C <- left_join(PAXINTEN_C, Mortality_2015_C, by = "SEQN") %>%
  left_join(Covariate_C, by=c("SEQN", "SDDSRVYR"))
AllAct_D <- left_join(PAXINTEN_D, Mortality_2015_D, by = "SEQN") %>%
  left_join(Covariate_D, by=c("SEQN", "SDDSRVYR"))

AllFlags_C <- left_join(Flags_C, Mortality_2015_C, by = "SEQN") %>%
  left_join(Covariate_C, by=c("SEQN", "SDDSRVYR"))
AllFlags_D <- left_join(Flags_D, Mortality_2015_D, by = "SEQN") %>%
  left_join(Covariate_D, by=c("SEQN", "SDDSRVYR"))

## clean up the workspace for memory purposes
rm(list=c(paste0(c("PAXINTEN_", "Covariate_","Mortality_2015_","Flags_"),rep(LETTERS[3:4],each=4))))

## combine data for the two waves
AllAct   <- rbind.data.frame(AllAct_C,AllAct_D)
AllFlags <- rbind.data.frame(AllFlags_C,AllFlags_D)


## clean up the workspace again
rm(list=c("AllAct_C","AllAct_D","AllFlags_C","AllFlags_D"))
##############################################################################
##                                                                          ##
##  Section 1b: create new variables/relevel factor variables for analyses  ##
##                                                                          ##
##############################################################################

## Code year 5 mortality, NAs for individuals with follow up less than 5 years and alive
AllAct$yr5_mort <- AllFlags$yr5_mort <- as.integer(ifelse(AllAct$permth_exm/12 <= 5 & AllAct$mortstat == 1, 1,
                                                          ifelse(AllAct$permth_exm/12 < 5 & AllAct$mortstat == 0, NA, 0)))
## Create Age in years using the age at examination (i.e. when participants wore the device)
AllAct$Age <- AllFlags$Age <- AllAct$RIDAGEEX/12

## Re-level comorbidities to assign refused/don't know as not having the condition
## Note that in practice this does not affect many individuals, but it is an 
## assumption we're making.
levels(AllAct$CHD) <- levels(AllFlags$CHD) <- list("No" = c("No","Refused","Don't know"),"Yes" = c("Yes"))
levels(AllAct$CHF) <- levels(AllFlags$CHF) <- list("No" = c("No","Refused","Don't know"),"Yes" = c("Yes"))
levels(AllAct$Stroke) <- levels(AllFlags$Stroke) <- list("No" = c("No","Refused","Don't know"),"Yes" = c("Yes"))
levels(AllAct$Cancer) <- levels(AllFlags$Cancer) <- list("No" = c("No","Refused","Don't know"),"Yes" = c("Yes"))
levels(AllAct$Diabetes) <- levels(AllFlags$Diabetes) <- list("No" = c("No","Borderline","Refused","Don't know"),"Yes" = c("Yes"))

## Re-level education to have 3 levels and categorize don't know/refused to be missing
levels(AllAct$EducationAdult) <- levels(AllFlags$EducationAdult) <- list("Less than high school"= c("Less than 9th grade","9-11th grade"),"High school" = c("High school grad/GED or equivalent"),"More than high school" = c("Some College or AA degree", "College graduate or above"))

## Re-level alcohol consumption to include a level for "missing"
levels(AllAct$DrinkStatus) <- levels(AllFlags$DrinkStatus) <- c(levels(AllAct$DrinkStatus),"Missing alcohol")
AllAct$DrinkStatus[is.na(AllAct$DrinkStatus)] <- AllFlags$DrinkStatus[is.na(AllAct$DrinkStatus)] <- "Missing alcohol"

## Re-order columns so that activity and wear/non-wear flags are the last 1440 columns of our two
## data matrices. This is a personal preference and is absolutely not necessary.
act_cols <- which(colnames(AllAct) %in% paste0("MIN",1:1440))
oth_cols <- which(!colnames(AllAct) %in% paste0("MIN",1:1440))
AllAct   <- AllAct[,c(oth_cols,act_cols)]
AllFlags <- AllFlags[,c(oth_cols,act_cols)]
rm(list=c("act_cols","oth_cols"))

###########################################################
##                                                       ##
##  Section 2: Calcualte common accelerometery features  ##
##                                                       ##
###########################################################

## Assign just the activity and wear/non-wear flag data to matrices.
## This makes computing the features faster but is technically required.
act_mat  <- as.matrix(AllAct[,paste0("MIN",1:1440)])
flag_mat <- as.matrix(AllFlags[,paste0("MIN",1:1440)])

## replace NAs with 0s
## As described in the manuscript, this only affects 501 minutes for 1 day, for one subject
act_mat[is.na(act_mat)]   <- 0
flag_mat[is.na(flag_mat)] <- 0

AllAct$TAC   <- AllFlags$TAC   <- rowSums(act_mat)
AllAct$TLAC  <- AllFlags$TLAC  <- rowSums(log(1+act_mat))
AllAct$WT    <- AllFlags$WT    <- rowSums(flag_mat)
AllAct$ST    <- AllFlags$ST    <- rowSums(act_mat < 100)
AllAct$MVPA  <- AllFlags$MVPA  <- rowSums(act_mat >= 2020)

## calculate fragmentation measures
bout_mat <- apply(act_mat >= 100, 1, function(x){
  mat <- rle(x)
  sed <- mat$lengths[which(mat$values == FALSE)]
  act <- mat$length[mat$values == TRUE]
  
  sed <- ifelse(length(sed) == 0, NA, mean(sed))
  act <- ifelse(length(act) == 0, NA, mean(act))
  c(sed,act)
})

AllAct$SBout <- AllFlags$SBout <- bout_mat[1,]
AllAct$ABout <- AllFlags$ABout <- bout_mat[2,]
AllAct$SATP  <- AllFlags$SATP  <- 1/AllAct$SBout
AllAct$ASTP  <- AllFlags$ASTP  <- 1/AllAct$ABout
rm(list=c("act_mat","flag_mat","bout_mat"))

###########################################
##                                       ##
##  Section 3: Apply exclusion criteria  ##
##                                       ##
###########################################

## make dataframe with one row per individual to create table 1.
## Remove columns associated with activity to avoid any confusion.
table_dat <- AllAct[!duplicated(AllAct$SEQN),-which(colnames(AllAct) %in% c(paste0("MIN",1:1440),
                                                                            "TAC","TLAC","WT","ST","MVPA","SBout","ABout","SATP","ASTP"))]

## subset based on our age inclusion/exclusion criteria
## note that individuals age 85 and over are coded as NA
table_dat <- subset(table_dat, !(Age < 50 | is.na(Age)))

## get the SEQN (id variable) associated with individuals with fewer than 3 days accelerometer 
## wear time with at least 10 hours OR had their data quality/device calibration flagged by NHANES
keep_inx       <- exclude_accel(AllAct, AllFlags)
Act_Analysis   <- AllAct[keep_inx,]
Flags_Analysis <- AllFlags[keep_inx,]
nms_rm         <- unique(c(Act_Analysis$SEQN[-which(Act_Analysis$SEQN %in% names(table(Act_Analysis$SEQN))[table(Act_Analysis$SEQN)>=3])],setdiff(AllAct$SEQN,Act_Analysis$SEQN)))
rm(list=c("keep_inx"))

## Additional inclusion/exclusion criteria.
## Aside from mortality or accelerometer weartime, the only missingness is in
## Education (6) and BMI (35).
criteria_vec <- c("(is.na(table_dat$BMI_cat))",         # missing BMI
                  "(is.na(table_dat$EducationAdult))",  # missing education
                  "(table_dat$SEQN %in% nms_rm)",       # too few "good" days of accel data
                  "((!table_dat$eligstat %in% 1) | is.na(table_dat$mortstat) | 
                  is.na(table_dat$permth_exm))") # missing mortality data, or accidental death

## create matrix of pairwise missing data based on our exclusion criterial
tab_miss <- matrix(NA, ncol=length(criteria_vec), nrow=length(criteria_vec))
for(i in seq_along(criteria_vec)){
  for(j in seq_along(criteria_vec)){
    eval(parse(text=paste0("miss_cur <- which(", criteria_vec[i], "&", criteria_vec[j],")")))
    tab_miss[i,j] <- length(miss_cur)
    rm(list=c("miss_cur"))
  }
}
rownames(tab_miss) <- colnames(tab_miss) <- c("BMI","Education","Bad Accel Data","Mortality")
rm(list=c("i","j"))
## view missing data pattern
tab_miss


## add in column indicating exclusion:
##   Exclude = 1 indicates an individual does not meet our inclusion criteria
##   Exclude = 0 indicates an individual does meet our inclusion criteria
eval(parse(text=paste0("table_dat$Exclude <- as.integer(", paste0(criteria_vec,collapse="|"), ")")))

## Create our dataset for analysis with one row per subject
## containing only those subjects who meet our inclusion criteria.
data_analysis  <- subset(table_dat, Exclude == 0)
data_analysis$mortstat <- ifelse((data_analysis$ucod_leading %in% "004" & data_analysis$mortstat ==1),0,data_analysis$mortstat)
## get adjusted survey weights using the reweight_accel function
data_analysis  <- reweight_accel(data_analysis)

## Get activity/flag data for only those included participants AND who have 3 good days of data.
## Since we've already removed the "bad" days from Act_Analysis and Act_Flags,
## we need only subset based on subject ID now.
Act_Analysis   <- subset(Act_Analysis, SEQN %in% data_analysis$SEQN)
Flags_Analysis <- subset(Flags_Analysis, SEQN %in% data_analysis$SEQN)

## calculate subject specific averages of the accelerometry features
## using only the "good" days of data
act_var_nms <- c("TAC","TLAC","WT","ST","MVPA","SATP","ASTP")
for(i in act_var_nms){
  data_analysis[[i]] <- vapply(data_analysis$SEQN, function(x) mean(Act_Analysis[[i]][Act_Analysis$SEQN==x],na.rm=TRUE), numeric(1))
}

## verify there's no missingness in the rest of our predictors of interest
vars_interest <- c("Age", "Gender", "Race", "EducationAdult", "SmokeCigs", "DrinkStatus", "BMI_cat",
                   "Diabetes","CHF",  "CHD", "Stroke",
                   "Cancer", "MobilityProblem",
                   "permth_exm")

## clean up the workspace
rm(list=c("AllAct","AllFlags","i","criteria_vec","nms_rm","tab_miss"))
gc()

###### data for EDA
data_analysis$time <- data_analysis$permth_exm/12
data_eda = data_analysis

# number of participants
nrow(data_analysis) 
# number of deaths
sum(data_analysis$mortstat==1) 
# person years of follow up time.
sum(data_analysis$time)