#!/bin/env Rscript

## Brent McPherson
## 20181126
## import the Cam-CAN behavioral test scores into a single data.frame
##

##library(xlsx) ## wtf?
library(readxl)
library(Amelia)
library(mice)
library(data.table)

## path to scored behavioral results
datadir <- "data"

## load all summary files in the latest release00# folder (slow)
files <- list.files(datadir, pattern = "summary.txt", all.files = TRUE, full.names = TRUE, recursive = TRUE)

## ... just import them 1 by 1 and merge on CCID dropping useless values (check w/ READMEs)

##
## just read them in and make variable names unique to data set...
##

## Benton Faces
dat01 <- read.table(files[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 657)
colnames(dat01)[2:length(colnames(dat01))] <- paste("mem_bfs", colnames(dat01)[2:length(colnames(dat01))], sep = "_")

## cardio measures
dat02 <- read.table(files[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 588)
colnames(dat02)[2:length(colnames(dat02))] <- paste("mtr_crd", colnames(dat02)[2:length(colnames(dat02))], sep = "_")

## Cattell
dat03 <- read.table(files[3], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 660)
colnames(dat03)[2:length(colnames(dat03))] <- paste("att_cat", colnames(dat03)[2:length(colnames(dat03))], sep = "_")

## Ekman
## read in individual tasks from single file
## make the vars b/w tasks unique before merging
dat04a <- read.table(files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 10, nrows = 664)
colnames(dat04a)[2:length(colnames(dat04a))] <- paste("ang", colnames(dat04a)[2:length(colnames(dat04a))], sep = "_")

dat04b <- read.table(files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 686, nrows = 664)
colnames(dat04b)[2:length(colnames(dat04b))] <- paste("dis", colnames(dat04b)[2:length(colnames(dat04b))], sep = "_")

dat04c <- read.table(files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1362, nrows = 664)
colnames(dat04c)[2:length(colnames(dat04c))] <- paste("fea", colnames(dat04c)[2:length(colnames(dat04c))], sep = "_")

dat04d <- read.table(files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2038, nrows = 664)
colnames(dat04d)[2:length(colnames(dat04d))] <- paste("hap", colnames(dat04d)[2:length(colnames(dat04d))], sep = "_")

dat04e <- read.table(files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 2714, nrows = 664)
colnames(dat04e)[2:length(colnames(dat04e))] <- paste("sad", colnames(dat04e)[2:length(colnames(dat04e))], sep = "_")

dat04f <- read.table(files[4], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 3390, nrows = 664)
colnames(dat04f)[2:length(colnames(dat04f))] <- paste("sur", colnames(dat04f)[2:length(colnames(dat04f))], sep = "_")

dat04 <- merge(dat04a[c(1, 3, 5, 6, 9:15)], dat04b[c(1, 3, 5, 6, 9:15)], by = "CCID", all = TRUE)
dat04 <- merge(dat04, dat04c[c(1, 3, 5, 6, 9:15)], by = "CCID", all = TRUE)
dat04 <- merge(dat04, dat04d[c(1, 3, 5, 6, 9:15)], by = "CCID", all = TRUE)
dat04 <- merge(dat04, dat04e[c(1, 3, 5, 6, 9:15)], by = "CCID", all = TRUE)
dat04 <- merge(dat04, dat04f[c(1, 3, 5, 6, 9:15)], by = "CCID", all = TRUE)

colnames(dat04)[2:length(colnames(dat04))] <- paste("emt_ekm", colnames(dat04)[2:length(colnames(dat04))], sep = "_")
               
## emotional memory
dat05 <- read.table(files[5], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 330)
colnames(dat05)[2:length(colnames(dat05))] <- paste("emt_emm", colnames(dat05)[2:length(colnames(dat05))], sep = "_")

## emotional regulation
dat06 <- read.table(files[6], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 314)
colnames(dat06)[2:length(colnames(dat06))] <- paste("emt_emr", colnames(dat06)[2:length(colnames(dat06))], sep = "_")

## famous faces
dat07a <- read.table(files[7], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 10, nrows = 660)
colnames(dat07a)[2:length(colnames(dat07a))] <- paste("fac", colnames(dat07a)[2:length(colnames(dat07a))], sep = "_")

dat07b <- read.table(files[7], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 682, nrows = 660)
colnames(dat07b)[2:length(colnames(dat07b))] <- paste("nam", colnames(dat07b)[2:length(colnames(dat07b))], sep = "_")

dat07 <- merge(dat07a[c(1, 5:12)], dat07b[c(1, 5:12)], by = "CCID", all = TRUE)
colnames(dat07)[2:length(colnames(dat07))] <- paste("mem_ffc", colnames(dat07)[2:length(colnames(dat07))], sep = "_")

## force matches
dat08 <- read.table(files[8], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 328)
colnames(dat08)[2:length(colnames(dat08))] <- paste("mtr_frc", colnames(dat08)[2:length(colnames(dat08))], sep = "_")

## Hotel
dat09 <- read.table(files[9], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 657)
colnames(dat09)[2:length(colnames(dat09))] <- paste("att_hot", colnames(dat09)[2:length(colnames(dat09))], sep = "_")

## MEG
#dat10 <- read.table(files[10], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 1000)
## empty? empty.

## motor learning
dat11 <- read.table(files[11], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 318)
colnames(dat11)[2:length(colnames(dat11))] <- paste("mtr_mtl", colnames(dat11)[2:length(colnames(dat11))], sep = "_")

## MRI 
dat12 <- read.table(files[12], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 657)
colnames(dat12)[2:length(colnames(dat12))] <- paste("mtr_mri", colnames(dat12)[2:length(colnames(dat12))], sep = "_")

## picture priming 
dat13 <- read.table(files[13], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 651)
colnames(dat13)[2:length(colnames(dat13))] <- paste("lng_pic", colnames(dat13)[2:length(colnames(dat13))], sep = "_")
colnames(dat13)[1] <- "CCID"

## proverbs
dat14 <- read.table(files[14], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 655)
colnames(dat14)[2:length(colnames(dat14))] <- paste("lng_pvb", colnames(dat14)[2:length(colnames(dat14))], sep = "_")

## choice reaction time
dat15 <- read.table(files[15], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 2501)
colnames(dat15)[2:length(colnames(dat15))] <- paste("att_crt", colnames(dat15)[2:length(colnames(dat15))], sep = "_")

## simple reaction time 
dat16 <- read.table(files[16], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 2550)
colnames(dat16)[2:length(colnames(dat16))] <- paste("att_srt", colnames(dat16)[2:length(colnames(dat16))], sep = "_")

## Synsem - sentence completion
dat17 <- read.table(files[17], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 634)
colnames(dat17)[2:length(colnames(dat17))] <- paste("lng_syn", colnames(dat17)[2:length(colnames(dat17))], sep = "_")
colnames(dat17)[1] <- "CCID"

## tip of the tongue
dat18 <- read.table(files[18], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 656)
colnames(dat18)[2:length(colnames(dat18))] <- paste("lng_tot", colnames(dat18)[2:length(colnames(dat18))], sep = "_")
colnames(dat18)[1] <- "CCID"

## VSTM color 
dat19 <- read.table(files[19], header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 8, nrows = 659)
colnames(dat19)[2:length(colnames(dat19))] <- paste("mem_vst", colnames(dat19)[2:length(colnames(dat19))], sep = "_")

##
## merge the datasets by subject ID
## 

## merge the first two data sets, grab variable names
#datRaw <- merge(dat01[c(1, 7, 8)], dat02[c(1, 9:17)], by = "CCID", all = TRUE)
datRaw <- dat01[c(1, 7, 8)] ## subscores
#datSum <- merge(dat01[c(1, 9)], dat02[c(1, 19:21)], by = "CCID", all = TRUE)
datReg <- dat02[c(1, 7, 8, 9:18)] ## physio

datRaw <- merge(datRaw, dat03[c(1, 4:7)], by = "CCID", all = TRUE) ## subscores
#datSum <- merge(datSum, dat03[c(1, 8)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat04[c(1, 6:11, 16:21, 26:31, 36:41, 46:51, 56:61)], by = "CCID") ## dummy codes
#datSum <- merge(datSum, dat04[c(1, 2:5, 12:15, 22:25, 32:35, 42:45, 52:55)], by = "CCID")

datRaw <- merge(datRaw, dat05[c(1, 11:65)], by = "CCID", all = TRUE) 
#datSum <- merge(datSum, dat05[c(1:10)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat06[c(1, 4, 5:7, 9:11, 13:15, 17:19)], by = "CCID", all = TRUE) ## conditions
#datSum <- merge(datSum, dat06[c(1, 4, 8, 12, 16, 20)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat07, by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat07, by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat08[c(1:17)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat08[c(1:17)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat09[c(1:3)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat09[c(1:3)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat11[c(1, 6:25)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat11[c(1:5)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat12[c(1, 4, 11, 13, 15, 16, 18)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat12[c(1, 4:8)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat13[c(1, 9:72)], by = "CCID", all = TRUE) ## actually not correlated 
#datSum <- merge(datSum, dat13[c(1:8)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat14[c(1, 7)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat14[c(1, 7)], by = "CCID", all = TRUE)

## multiple redundant
## inv = inverted to normalize better; trim3 = removed above 3SD outliers; _## finger used
d15.indx <- colnames(dat15) %like% "inv" & colnames(dat15) %like% "trim3" 
d15.drop <- colnames(dat15) %like% "cv" | colnames(dat15) %like% "median" | colnames(dat15) %like% "all"
##colnames(dat15)[d15.indx & !d15.drop]

##datRaw <- merge(datRaw, dat15[c(1, 2, 15:132)], by = "CCID", all = TRUE)
datRaw <- merge(datRaw, dat15[c(1, 6, 7, 8, 10, which(d15.indx & !d15.drop))], by = "CCID", all = TRUE)
## need to fix global merge if I use it at all
##datSum <- merge(datSum, dat15[c(1, 2, 6:14)], by = "CCID", all = TRUE)

##datRaw <- merge(datRaw, dat16[c(1, 2, 13:26)], by = "CCID", all = TRUE)
d16.indx <- colnames(dat16) %like% "inv" & colnames(dat16) %like% "trim3"
d16.drop <- colnames(dat16) %like% "cv" | colnames(dat16) %like% "median"

datRaw <- merge(datRaw, dat16[c(1, 5:7, which(d16.indx & !d16.drop))], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat16[c(1, 2, 4:12)], by = "CCID", all = TRUE)

## keep specific responses to all individual conditions, not overall
d17.keep <- colnames(dat17) %like% "pYes" | colnames(dat17) %like% "reITRTvalid" |
    colnames(dat17) %like% "reITRTyes" | colnames(dat17) %like% "reITRTno" |
    colnames(dat17) %like% "pValid" | colnames(dat17) %like% "pAnticipation"
d17.drop <- colnames(dat17) %like% "overall"

datRaw <- merge(datRaw, dat17[c(1, 2, which(d17.keep & !d17.drop))], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat17[c(1:9)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat18[c(1:7)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat18[c(1, 8)], by = "CCID", all = TRUE)

datRaw <- merge(datRaw, dat19[c(1:13)], by = "CCID", all = TRUE)
#datSum <- merge(datSum, dat19[c(1, 14, 15)], by = "CCID", all = TRUE)

##
## import participant_data.tsv
##

## read in data
prtd <- read.table("participant_data.tsv", header = TRUE, stringsAsFactors = FALSE)

## fix variable names / coding
colnames(prtd)[1] <- "CCID"
colnames(prtd)[5] <- "male"
prtd$male <- ifelse(prtd$male == 1, 1, 0)

## subset out useful variables for analysis
prtd <- prtd[, c(1:3, 5, 6)]

## merge together the other nuisance regressors to merge to the front of the data set
datReg <- merge(prtd, datReg, by = "CCID", all = TRUE)

##
## load in home interview questions
##

#hint <- read.xlsx("approved_data.xlsx", na.strings = "NaN")
hint <- read_xlsx("approved_data.xlsx", 1, na = "NaN")
hint <- data.frame(hint)

## code gender sanely
colnames(hint)[2] <- "male"
hint$male <- ifelse(hint$male == "F", 0, 1)
hint$homeint_v45 <- ifelse(hint$homeint_v45 == 1, 1, 0)

## png(filename='figs/missmap.png')
## missmap(hint)
## dev.off()

## filter is some proportion are missing, regardless of what they are?

## index 4     : mmse total 
## index 5     : handedness
## index 7     : number of correct 7's (basic cognitive test)
## index 9     : marital status
## index 10    : home status
## index 13    : funds accomadations
## index 15    : owned / rented
## index 16    : how many live here
## index 18    : average income of residents
## index 24    : bilingual from birth
## index 32    : ever had paid work?
## index 37-40 : manage multiple employees
## index 42:44 : type of work
## index 45:48 : parters sex / type of work
## index 64:65 : education completed / when
## index 77:79 : have children
## index 80:89 : family / friend interactions

## homeint_v241-261 : smoke?
## homeint_v262-289 : diet
## homeint_v290-299 : alcohol
## homeint_v300-320 : drug abuse
## homeint_v321-347 : overall health, vision, hearing
## homeint_v339+ : health condition / age diagnosed
## homeint_v481+ : sleep
## homeint_v504+ : handedness

## homeint_v515+ : story recall
## homeint_v528+ : anxiety / affect
## homeint_v542+ : fallen
## homeint_v565+ : self-care
## homeint_v~600 : medical history

## dat20 <- hint[c(1, 4, 5, 7, 9, 10, 13, 15, 16, 18, 24, 32, 37:40, 42:44, 45:48, 64, 65, 77:79, 80:89)]

##
## home interview questions to keep
##

## paste together question stem and numbers to find variable names to keep as nuissance regressors
## drop 15 (income) b/c they don't translate the letters to numbers in the key
## drop 73 (qualifications) b/c there's no easy way to combine multiple responses
## drop 96 (meetings) b/c there's no easy way to combine multiple responses + frequency added
hq <- paste('homeint_v', c(5, 6, 10:13, 19, 21, 22, 29, 34:48, 74, 86:95, 97, 98,
                          230:240, 321:323, 325, 349, 354, 355, 461, 481:486, 492, 493,
                          496:503, 599, 600:602, 604, 606, 607, 609, 611, 612, 614:616), sep = "")
hm <- colnames(hint)[colnames(hint) %like% "mmse"]
hq <- c(hq, paste('homeint_', c('recheck', 'sevens_n'), sep = ""), hm)
hq <- c(c("CCID", 'male'), hq, hm)

## home interview questionaire
hqi <- colnames(hint) %in% hq
hqd <- hint[, hqi]

## set "." values to missing - not sure how this time was converted...
hqd$homeint_v481 <- as.numeric(ifelse(hqd$homeint_v481 == ".", NA, hqd$homeint_v481))
hqd$homeint_v483 <- as.numeric(ifelse(hqd$homeint_v483 == ".", NA, hqd$homeint_v483))
hqd$homeint_v484 <- as.numeric(ifelse(hqd$homeint_v484 == ".", NA, hqd$homeint_v484))

## rename into usefull variables
colnames(hqd) <- gsub("homeint", "hint", colnames(hqd))
colnames(hqd)[3] <- "hint_mmse_interview"
colnames(hqd)[4] <- "hint_married"
colnames(hqd)[5] <- "hint_accommodation"
colnames(hqd)[6] <- "hint_accom_fund"
colnames(hqd)[7] <- "hint_accom_years"
colnames(hqd)[8] <- "hint_accom_own"
colnames(hqd)[9] <- "hint_accom_nres"
colnames(hqd)[10] <- "hint_lang_1"
colnames(hqd)[11] <- "hint_bilingual"
colnames(hqd)[12] <- "hint_lang_2"
colnames(hqd)[13] <- "hint_paid_work"
colnames(hqd)[14] <- "hint_age_retired"
colnames(hqd)[15] <- "hint_self_employed"
colnames(hqd)[16] <- "hint_n_employees"
colnames(hqd)[17] <- "hint_supervisor"
colnames(hqd)[18] <- "hint_n_years_worked"
colnames(hqd)[19] <- "hint_hours_per_week"
colnames(hqd)[20] <- "hint_shift_work"
colnames(hqd)[21] <- "hint_night_shift"
colnames(hqd)[22] <- "hint_partner_male"
colnames(hqd)[23] <- "hint_partner_employ"
colnames(hqd)[24] <- "hint_partner_paid"
colnames(hqd)[25] <- "hint_partner_age_retired"
colnames(hqd)[26] <- "hint_age_edu_finish"
colnames(hqd)[27] <- "hint_have_kids"
colnames(hqd)[28] <- "hint_kids_living"
colnames(hqd)[29] <- "hint_kids_dead"
colnames(hqd)[30] <- "hint_relatives_nearby"
colnames(hqd)[31] <- "hint_relatives_see_freq"
colnames(hqd)[32] <- "hint_relatives_speak_freq"
colnames(hqd)[33] <- "hint_relatives_txt_email"
colnames(hqd)[34] <- "hint_visit_freq"
colnames(hqd)[35] <- "hint_friends_speak"
colnames(hqd)[36] <- "hint_friends_txt_email"
colnames(hqd)[37] <- "hint_friends_community"
colnames(hqd)[38] <- "hint_see_neighbors"
colnames(hqd)[39] <- "hint_recall_mem_prob"
colnames(hqd)[40] <- "hint_recall_what_day"
colnames(hqd)[41] <- "hint_recall_unkreps"
colnames(hqd)[42] <- "hint_recall_unkdeath"
colnames(hqd)[43] <- "hint_recall_unk_mnth_yr"
colnames(hqd)[44] <- "hint_recall_unkaction"
colnames(hqd)[45] <- "hint_recall_getlost"
colnames(hqd)[46] <- "hint_recall_cantfind"
colnames(hqd)[47] <- "hint_recall_read"
colnames(hqd)[48] <- "hint_recall_readtv"
colnames(hqd)[49] <- "hint_recall_impact"
colnames(hqd)[50] <- "hint_health_est"
colnames(hqd)[51] <- "hint_health_last_year"
colnames(hqd)[52] <- "hint_health_chronic_iss"
colnames(hqd)[53] <- "hint_health_wght_chng"
colnames(hqd)[54] <- "hint_health_high_bp"
colnames(hqd)[55] <- "hint_health_bp_untrt"
colnames(hqd)[56] <- "hint_health_choles"
colnames(hqd)[57] <- "hint_health_bone_frac"
colnames(hqd)[58] <- "hint_sleep_bed_time"
colnames(hqd)[59] <- "hint_sleep_fall_time"
colnames(hqd)[60] <- "hint_rise_time"
colnames(hqd)[61] <- "hint_sleep_hours"
colnames(hqd)[62] <- "hint_sleep_win_30min"
colnames(hqd)[63] <- "hint_sleep_interrupt"
colnames(hqd)[64] <- "hint_bad_dreams"
colnames(hqd)[65] <- "hint_sleep_pain"
colnames(hqd)[66] <- "hint_sleep_quality"
colnames(hqd)[67] <- "hint_sleep_meds"
colnames(hqd)[68] <- "hint_diff_stay_awake"
colnames(hqd)[69] <- "hint_no_enthusiasm"
colnames(hqd)[70] <- "hint_bed_partner"
colnames(hqd)[71] <- "hint_tired_morning"
colnames(hqd)[72] <- "hint_tired_afternoon"
colnames(hqd)[73] <- "hint_tired_evening"
colnames(hqd)[74] <- "hint_adopted"
colnames(hqd)[75] <- "hint_twin_triplet"
colnames(hqd)[76] <- "hint_dad_alive"
colnames(hqd)[77] <- "hint_dad_age"
colnames(hqd)[78] <- "hint_dad_age_died"
colnames(hqd)[79] <- "hint_mom_alive"
colnames(hqd)[80] <- "hint_mom_age"
colnames(hqd)[81] <- "hint_mom_age_died"
colnames(hqd)[82] <- "hint_nbrothers"
colnames(hqd)[83] <- "hint_nsisters"
colnames(hqd)[84] <- "hint_health_meds"
colnames(hqd)[85] <- "hint_health_prescribed"
colnames(hqd)[86] <- "hint_health_otc"

## drop variables unlikely to be useful after comparing w/ other data sets
hqd <- hqd[, c(-2, -3, -10, -12, -(39:49), -56, -(58:73), -84, -87, -88, -117)]
## drop variables that overlap w/ add summaries or are exlusion critieria

## split into actual values for estimation + regressors

## pull values for inclusion in CCA, recode
hqe <- hqd[, c(1, 81, 55:80, 8, 22, 24:38, 42:53)]
colnames(hqe)[2:28] <- gsub("hint", "cogs", colnames(hqe)[2:28])
colnames(hqe)[33:42] <- gsub("hint", "soc", colnames(hqe)[33:42])

## pull variables to regress out
hqr <- hqd[, c(1:7, 9:21, 39:41)]
colnames(hqr)[2:dim(hqr)[2]] <- paste("reg", colnames(hqr)[2:dim(hqr)[2]], sep = "_")

##
## physical activity self report
##

## exercise / activity measures
eqi <- colnames(hint) %like% 'epaq'
eqi[1] <- TRUE
eqd <- hint[, eqi] ## activity / mobility measures?
colnames(eqd)[2:dim(eqd)[2]] <- paste("reg", colnames(eqd)[2:dim(eqd)[2]], sep = "_")

##
## life experiences self report
##

## scq_*  : life experiences / employment
sqi <- colnames(hint) %like% 'scq_'
sqi[1] <- TRUE
sqd <- hint[, sqi]

## only keep the numerically coded variables (these aren't worth parsing from strings to numbers)
sqk <- (sapply(sqd, class) == "numeric") | (sapply(sqd, class) == "logical")
names(sqk) <- NULL
sqk[1] <- TRUE

## create the final subset of variables
leq <- sqd[, sqk]

## change the variable label to life experiences questions
colnames(leq) <- gsub("scq", "reg_leq", colnames(leq))

##
## composite scores (other cognitive scales)
##

## additional_* : derived question scores from raw questions
adi <- colnames(hint) %like% 'additional'
adi[1] <- TRUE
add <- hint[, adi]

add$additional_student_iv <- ifelse(add$additional_student_iv == "Yes", 1, 0)

add$additional_drug_severity <- ifelse(add$additional_drug_severity == "Low", 1,
                                       ifelse(add$additional_drug_severity == "Intermediate", 2, 3))

add$additional_drug_intervention <- ifelse(add$additional_drug_intervention == "Brief Intervention", 1,
                                    ifelse(add$additional_drug_intervention == "Outpatient (Intensive)", 2, 3))

add$additional_smoker <- ifelse(add$additional_smoker == "Never smoked", 1,
                         ifelse(add$additional_smoker == "Past smoker", 2,
                         ifelse(add$additional_smoker == "Non current, less than 100", 3,
                         ifelse(add$additional_smoker == "Non current, at least 100", 4,
                         ifelse(add$additional_smoker == "Occasional current, less than 100", 5,
                         ifelse(add$additional_smoker == "Occasional current, at least 100", 6, 7))))))

add$additional_alcohol <- ifelse(add$additional_alcohol == "Non drinker", 1,
                          ifelse(add$additional_alcohol == "Past drinker", 2,
                          ifelse(add$additional_alcohol == "Occasional drinker", 3,
                          ifelse(add$additional_alcohol == "Drinks one/three time monthly", 4,
                          ifelse(add$additional_alcohol == "Drinks once/twice weekly", 5,
                          ifelse(add$additional_alcohol == "Drinks three/four times weekly", 6, 7))))))

add$additional_HADS_anx_category <- ifelse(add$additional_HADS_anx_category == "Normal", 1,
                                    ifelse(add$additional_HADS_anx_category == "Mild", 2,
                                    ifelse(add$additional_HADS_anx_category == "Moderate", 3, 4)))

add$additional_HADS_dep_category <- ifelse(add$additional_HADS_dep_category == "Normal", 1,
                                    ifelse(add$additional_HADS_dep_category == "Mild", 2,
                                    ifelse(add$additional_HADS_dep_category == "Moderate", 3, 4)))

add$additional_qual <- ifelse(add$additional_qual == "None", 1,
                       ifelse(add$additional_qual == "GCSE/O-level", 2,
                       ifelse(add$additional_qual == "A-level", 3, 4)))
                              
add$additional_SC_r <- ifelse(add$additional_SC_r == "-", 1,
                       ifelse(add$additional_SC_r == "I", 2,
                       ifelse(add$additional_SC_r == "II", 3,
                       ifelse(add$additional_SC_r == "IIIM", 4,
                       ifelse(add$additional_SC_r == "IIIN", 5,
                       ifelse(add$additional_SC_r == "IV", 6, 7))))))
                              
add$additional_SC_pt <- ifelse(add$additional_SC_pt == "-", 1,
                        ifelse(add$additional_SC_pt == "I", 2,
                        ifelse(add$additional_SC_pt == "II", 3,
                        ifelse(add$additional_SC_pt == "IIIM", 4,
                        ifelse(add$additional_SC_pt == "IIIN", 5,
                        ifelse(add$additional_SC_pt == "IV", 6, 7))))))

add$additional_SEG_r <- as.numeric(add$additional_SEG_r)

## shorten variable names
colnames(add) <- gsub("additional", "comp", colnames(add))

colnames(add)[17] <- "comp_acer_attention"
colnames(add)[18] <- "comp_acer_memory"
colnames(add)[19] <- "comp_acer_fluencies"
colnames(add)[20] <- "comp_acer_language"
colnames(add)[21] <- "comp_acer_visuospatial"

colnames(add)[23] <- "comp_psqi_sleep_quality"
colnames(add)[24] <- "comp_psqi_sleep_latency"
colnames(add)[25] <- "comp_psqi_sleep_duration"
colnames(add)[26] <- "comp_psqi_sleep_efficiency"
colnames(add)[27] <- "comp_psqi_sleep_disturbance"
colnames(add)[28] <- "comp_psqi_sleep_medication"
colnames(add)[29] <- "comp_psqi_sleep_day_disturb"

## drop student IV (nearly singular)
add <- add[, -2]

##
## merge for export of behavioral data
##

## merge EPAQ / LEQ
demReg <- merge(eqd, leq, by = "CCID", all = TRUE)

## merge in demographic regressors
demReg <- merge(demReg, hqr, by = "CCID", all = TRUE)

## merge in other regressors
outReg <- merge(datReg, demReg, by = "CCID", all = TRUE)

## separate component scores from composites in add for HADS / ACER / PSQI
##adr <- add[, c(1:8, 9, 10, 13:20, 22:31, 33:36)] ## raw scores - too many dropped
##ads <- add[, c(1:8, 11:15, 21, 32, 33:36)] ## summary scores - too many excluded
ado <- add[, c(1:10, 13:15, 21:36)] ## just use these

## separate component scores from composites in hint for MMSE
reg <- merge(ado, hqe, by = "CCID", all = TRUE)
## reg <-reg[, c(1, 33:59, 2:32, 60:63, 74:88, 64:73)] ## indices sorting fro adr
reg <-reg[, c(1, 30, 2:29, 57:60, 71:85, 61:70)] 
    
## merge demographics / screening w/ behavior (datRaw)
outRaw <- merge(reg, datRaw, by = "CCID", all = TRUE)

##
## sort and merge final data sets
##

## load the subject ID's w/ scan data
## ccid <- read.table("ccid_old.txt", stringsAsFactors = FALSE)
## colnames(ccid) <- "CCID"

## merge to fill in missing
outAll <- merge(outReg, outRaw, by = "CCID", all = TRUE)
## outAll <- outAll[outAll$CCID %in% ccid$CCID, ]

## subset down to the subjects that have networks
outReg <- outAll[, 1:544]
outRaw <- outAll[, c(1, 545:938)]

## convert logical values in regressors to numeric values
regCols <- sapply(outReg, is.logical)
outReg[, regCols] <- lapply(outReg[, regCols], as.numeric)

##
## save .csv files of data
##

## write the data
write.table(outReg[, -1], 'camcan_reg_data.csv', quote = FALSE, sep = ",", na = "NaN",
            row.names = FALSE, col.names = FALSE)
write.table(outRaw[, -1], 'camcan_raw_data.csv', quote = FALSE, sep = ",", na = "NaN",
            row.names = FALSE, col.names = FALSE)

## write the subject IDs
write.table(outReg[, 1], 'camcan_reg_id.csv', quote = FALSE, sep = ",", na = "NaN",
            row.names = FALSE, col.names = FALSE)
write.table(outRaw[, 1], 'camcan_raw_id.csv', quote = FALSE, sep = ",", na = "NaN",
            row.names = FALSE, col.names = FALSE)

## write the variable names
write.table(colnames(outReg)[-1], 'camcan_reg_var_names.csv', quote = FALSE, sep = ",", na = "NaN",
            row.names = FALSE, col.names = FALSE)
write.table(colnames(outRaw)[-1], 'camcan_raw_var_names.csv', quote = FALSE, sep = ",", na = "NaN",
            row.names = FALSE, col.names = FALSE)

##
## Original: load in subject ID's w/ scan data
##

## ## load the subject ID's w/ scan data
## ccid <- read.table("ccid.txt", stringsAsFactors = FALSE)
## colnames(ccid) <- "CCID"

## ## merge the nuisance regressors for removal from the front of the data
## datRaw <- merge(datReg, datRaw, by = "CCID", all = TRUE)
## #datSum <- merge(datReg, datSum, by = "CCID", all = TRUE)

## ## create output data that is all behavior intersecting w/ scans
## outRaw <- datRaw[datRaw$CCID %in% ccid$CCID, ]
## #outSum <- datSum[datSum$CCID %in% ccid$CCID, ]
## # tmpReg <- datReg[datReg$CCID %in% ccid$CCID, ]

## # tmp <- merge(tmpReg, datSum, by = "CCID", all = TRUE)
## # outReg <- tmp[tmp$CCID %in% ccid$CCID, 1:4]

## missmap(outRaw)
## #missmap(outSum)
## # missmap(outReg)

## outRawID <- outRaw$CCID
## outSumID <- outSum$CCID
## #outRegID <- outReg$CCID

## ## SET MISSING AS NAN

## ## write out subject IDs
## write.table(outRawID, 'camcan_594_raw_id.csv', quote = FALSE, sep = ",", na = "NaN",
##             row.names = FALSE, col.names = FALSE)
## write.table(outSumID, 'camcan_594_sum_id.csv', quote = FALSE, sep = ",", na = "NaN",
##             row.names = FALSE, col.names = FALSE)
## #write.table(outRegID, 'camcan_594_reg_id.csv', quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

## ## write out variable names
## write.table(colnames(outRaw)[2:length(colnames(datRaw))], 'camcan_594_raw_vars.csv', quote = FALSE, sep = ",",
##                      na = "NaN", row.names = FALSE, col.names = FALSE)
## write.table(colnames(outSum)[2:length(colnames(datSum))], 'camcan_594_sum_vars.csv', quote = FALSE, sep = ",",
##                      na = "NaN"row.names = FALSE, col.names = FALSE)

## ## write out data
## write.table(outRaw[, -1], 'camcan_594_raw.csv', quote = FALSE, sep = ",", na = "NaN",
##             row.names = FALSE, col.names = FALSE)
## write.table(outSum[, -1], 'camcan_594_sum.csv', quote = FALSE, sep = ",", na = "NaN",
##             row.names = FALSE, col.names = FALSE)
## #write.table(outReg[, -1], 'camcan_594_reg.csv', quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)

## ## pull age from data before processing once ID's are meged w/ the networks





##
## test mice MI
##

## only use variables w/ a minimum correlation to impute values
#pred <- quickpred(outRaw, 0.3) ## failed to estimate

## initialize
#ini <- mice(outRaw, method = "cart", maxit = 0, seed = 66045)

## run a single imputation for testing
tmp <- mice(outRaw, m = 1, method = "cart", maxit = 5, seed = 66045, print = FALSE)

## create the completed data set
out <- complete(tmp)

## save imputed complete data set to disk
write.table(out[, -1], 'camcan_594_raw_imputed.csv', quote = FALSE, sep = ",",
            row.names = FALSE, col.names = FALSE)

## see if Amelia returns a single data set to use
?amelia

?mi.meld ## I think this is merging results, not data sets
## https://rdrr.io/cran/Amelia/man/mi.meld.html

tmp <- amelia(outRaw[, -1]) ## fails w/ multiple vars...

