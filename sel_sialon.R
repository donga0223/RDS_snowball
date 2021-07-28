library(tidyverse)
library(haven)
library(expss)
library(RDS)
library(igraph)
#sialon<- read_dta(file = "/Users/dongahkim/Dropbox/RDS_snowball/Sialon II_RDS_KG.dta") 
sialon<- read_dta(file = "/Users/dongahkim/OneDrive/RDS_snowball/Data/RDS_final_Krista.dta", encoding = "latin1") 
sialon$rec.id <- substr(sialon$RDS_CODE,1,nchar(sialon$RDS_CODE)-1)
sialon$center <- substr(sialon$RDS_CODE,1,2)

sialon$q16[which(is.na(sialon$q16))] <- 0
sialon$q18[which(is.na(sialon$q18))] <- 0
sialon$degree <- sialon$SN_FINAL

sialon$rec.id[which(sialon$rec.id=="IT")] <- 0
sialon$rec.id[which(sialon$rec.id=="SK")] <- 0
sialon$rec.id[which(sialon$rec.id=="RO")] <- 0
sialon$rec.id[which(sialon$rec.id=="LT")] <- 0

q10 <- ifelse(sialon$q10==1,1,0)
q10[is.na(sialon$q10)] <- 0


q11a <- ifelse(sialon$q11a==1,1,0)
q11a[is.na(sialon$q11a)] <- 0
q11b<- ifelse(sialon$q11b==1,1,0)
q11b[is.na(sialon$q11b)] <- 0
q11c <- ifelse(sialon$q11c==1,1,0)

q11c[is.na(sialon$q11c)] <- 0
q11d <- ifelse(sialon$q11d==1,1,0)
q11d[is.na(sialon$q11d)] <- 0

q14 <- ifelse(sialon$q14==1,1,0)
q14[is.na(sialon$q14)] <- 0

q33 <- ifelse(sialon$q33==1,1,0)
q33[is.na(sialon$q33)] <- 0

q12 <- ifelse(sialon$q12==1,1,0)
q12[is.na(sialon$q12)] <- 0

q35 <- ifelse(sialon$q35==1,1,0)
q35[is.na(sialon$q35)] <- 0


homosexual <- ifelse(sialon$q45==1,1,0)
sel.sialon <- data.frame(RDS_CODE = sialon$RDS_CODE, age = sialon$age, HIV = sialon$HIV
                         , current_live =  sialon$q3, where_live = sialon$q4
                         , HSU = q10, q11a, q11b, q11c, q11d, rec.id = sialon$rec.id
                         , degree = sialon$degree, center_region = sialon$center
                         , HIV_testing_history = q14, other_testing_history = q12
                         , ART_coverage = q33, howmany_nonsteady_malepartners = sialon$q16
                         , howmany_nonsteady_malepartners_unprotected = sialon$q18
                         , injected_drug = q35, current_occupation = sialon$q41
                         , yourself = sialon$q45, homosexual
                         , COUPON1 = sialon$COUPON1, COUPON2 = sialon$COUPON2
                         , COUPON3 = sialon$COUPON3)

write.csv(sel.sialon, "sel.sialon.csv")

#write.csv(sel.sialon, "sel.sialon.csv")
