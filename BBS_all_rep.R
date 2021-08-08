library(tidyverse)
library(haven)
library(expss)
library(RDS)
library(igraph)
library(base)


df <- read_sav("/Users/dongahkim/OneDrive/RDS_snowball/Data/462.GF_BBS_PWID_2016_17_7cities.sav")
dim(df)
#names(df)
#df[["I4_1"]]
rec.id <-  substr(df$CouponN,1,nchar(df$CouponN)-1)
rec.id[which(rec.id == "")] <- 0
rec.id <- as.numeric(rec.id)
received_syringes_12 <- ifelse(df$I4_1 == 1, 1, 0)  ## 0 은?? 1 yes, 2 no
gender <- ifelse(df$A6_gender == 1, 1, 0)
df1 <- data.frame(ID = df$ID, city = df$city, degree = df$S_7, RDS_CODE = df$CouponN
                  , rec.id
                  #, Coupon1 = df$Coupon1, Coupon2 = df$Coupon2, Coupon3 = df$Coupon3
                  #, C1, C2, C3
                  , age = df$A5_age, gender, HIV = df$HIV_1, HCV = df$Hep_C
                  , tested_know = df$Tested_know_2, sage_inj_last = df$safe_inj_last_11
                  , cond_use_last_sex = df$cond_use_last_sex_12_rec
                  , received_syringes_12)

cityname <- c("Tbilisi", "Gori", "Telavi", "Zugdidi", "Batumi", "Kutaisi", "Rustavi")


get.RDS.data <- function(data){
  RDS.data <- as.rds.data.frame(data, id = "RDS_CODE"
                                , recruiter.id = "rec.id"
                                , network.size = "degree")
  
  RDS.data$C3 <- RDS.data$C2 <- RDS.data$C1 <- NA
  
  for(i in 1:dim(RDS.data)[1]){
    if(length(which(RDS.data$rec.id==RDS.data$RDS_CODE[i]))==0){
      RDS.data[i,"C1"] <- RDS.data[i,"C2"] <- RDS.data[i,"C3"] <- NA
    }else if(length(which(RDS.data$rec.id==RDS.data$RDS_CODE[i]))==1){
      RDS.data[i,"C1"] <- RDS.data$RDS_CODE[which(RDS.data$rec.id==RDS.data$RDS_CODE[i])]
      RDS.data[i,"C2"] <- RDS.data[i,"C3"] <- NA
    }else if(length(which(RDS.data$rec.id==RDS.data$RDS_CODE[i]))==2){
      RDS.data[i,"C1"] <- RDS.data$RDS_CODE[which(RDS.data$rec.id==RDS.data$RDS_CODE[i])][1]
      RDS.data[i,"C2"] <- RDS.data$RDS_CODE[which(RDS.data$rec.id==RDS.data$RDS_CODE[i])][2]
      RDS.data[i,"C3"] <- NA
    }else if(length(which(RDS.data$rec.id==RDS.data$RDS_CODE[i]))==3){
      RDS.data[i,"C1"] <- RDS.data$RDS_CODE[which(RDS.data$rec.id==RDS.data$RDS_CODE[i])][1]
      RDS.data[i,"C2"] <- RDS.data$RDS_CODE[which(RDS.data$rec.id==RDS.data$RDS_CODE[i])][2]
      RDS.data[i,"C3"] <- RDS.data$RDS_CODE[which(RDS.data$rec.id==RDS.data$RDS_CODE[i])][3]
    }
  }
  return(RDS.data)
}

pdf("BBS.RDSplot.2016_17_7cities.pdf")
for(i in 1:length(cityname)){
  assign(paste("df1.", cityname[i], sep=""), df1[which(df1$city == i),])
  assign(paste("RDS.", cityname[i], sep=""),  get.RDS.data(get(paste("df1.", cityname[i], sep=""))))
  
  reingold.tilford.plot(get(paste("RDS.", cityname[i], sep="")), vertex.color = "HIV", vertex.size = "received_syringes_12"
                        , vertex.size.range = c(2,5), main = cityname[i])

}
dev.off()


source("/Users/dongahkim/OneDrive/RDS_snowball/Rcode/snowball_RDS_resampling_ftn.R")
BBS.est <- function(data, RDS.data){
  res <- data.frame()
  aa <- data %>% 
    summarise(n = n(), mean(received_syringes_12, na.rm = TRUE), mean(age, na.rm = TRUE)
              , sum(gender[which(gender == 1)], na.rm = TRUE)/(n())
              , mean(HIV, na.rm = TRUE), mean(HCV, na.rm = TRUE)
              , mean(tested_know, na.rm = TRUE)
              , mean(sage_inj_last, na.rm = TRUE)
              , mean(cond_use_last_sex, na.rm = TRUE))
  
  
  
  res <- as.data.frame(aa)
  names(res) <- c("n", "received_syringes_12", "age", "gender"
                  , "HIV", "HCV", "tested_know", "sage_inj_last"
                  , "cond_use_last_sex")
                  
  
  myvar <- c("received_syringes_12", "gender"
             , "HIV", "HCV", "tested_know", "sage_inj_last"
             , "cond_use_last_sex")
  for(k in 1:length(myvar)){
    assign(myvar[k], RDS.II.estimates(rds.data=RDS.data ,outcome.variable=myvar[k]))
  }
  
  res.RDS <- c()
  for(k in 1:length(myvar)){
    res.RDS <-  as.data.frame(cbind(res.RDS, get(myvar[k])$interval[1,1:3]))
  }
  
  colnames(res.RDS) <- myvar
  #colnames(res.RDS) <- c("Estimate", "95% LCI", "95% UCI")
  
  
  return(list(res = res, res.RDS = res.RDS))
}


#snow.data <- snowball.wave1(RDS.data = RDS.Tbilisi, sample.size = dim(RDS.data)[1]
#                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "received_syringes_12")
#aa <- RDS.sample.dongah(RDS.data = RDS.Tbilisi, sample.size = dim(RDS.data)[1]
#                              , number.of.coupons = 3, sample.with.replacement = TRUE)
  


snowball.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway, seed.var = seed.var){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data, seed.var = seed.var)
  snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")
  
  HSU <- data.frame()
  aa <- BBS.est(snow.data, snow.RDS)
  HSU <- as.data.frame(cbind(aa$res, aa$res.RDS[1,]))
  
  for(i in 2:iter){
    set.seed(setseed+i)
    snow.data <- snowball.sampleway(RDS.data, seed.var = seed.var)
    snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
    aa <- BBS.est(snow.data, snow.RDS)
    HSU <- rbind(HSU, as.data.frame(cbind(aa$res, aa$res.RDS[1,])))
    print(i)
  }
  return(HSU) 
}

for(i in 1:length(cityname)){
  assign(paste(cityname[i], ".wave2.rep", sep=""), snowball.rep(RDS.data = get(paste("RDS.", cityname[i], sep="")), 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "received_syringes_12"))
  print(i)
  save.image("BBS_all_rep.RData")
}



RDS.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway = RDS.sample.dongah){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data)
  snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")
  
  HSU <- data.frame()
  aa <- BBS.est(snow.data, snow.RDS)
  HSU <- as.data.frame(cbind(aa$res, aa$res.RDS[1,]))
  
  
  for(p in 2:iter){
    set.seed(setseed+p)
    snow.data <- snowball.sampleway(RDS.data)
    snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
    aa <- BBS.est(snow.data, snow.RDS)
    HSU <- rbind(HSU, as.data.frame(cbind(aa$res, aa$res.RDS[1,])))
    print(p)
  }
  return(HSU) 
}

for(i in 1:length(cityname)){
  assign(paste(cityname[i], ".RDS.rep", sep=""), RDS.rep(RDS.data = get(paste("RDS.", cityname[i], sep="")), 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah))
  print(i)
  save.image("BBS_all_rep.RData")
}