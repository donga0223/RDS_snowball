library(tidyverse)
library(haven)
library(expss)
library(RDS)
library(igraph)
library(base)


#source("SPRTBA_function.R")
sel.sialon <- read.csv("sel.sialon.csv")
#sel.sialon <- read.csv("/Users/dongahkim/Dropbox/RDS_snowball/sel.sialon.csv")

sialon.IT <- sel.sialon[which(sel.sialon$center == "IT"),]
sialon.SK <- sel.sialon[which(sel.sialon$center == "SK"),]
sialon.RO <- sel.sialon[which(sel.sialon$center == "RO"),]
sialon.LT <- sel.sialon[which(sel.sialon$center == "LT"),]

print("1.done")

data <- sialon.IT
data[which(data$degree> dim(data)[1]),"degree"] <- dim(data)[1]/2

print("1.done")

########### as.RDS.data
get.RDS.data <- function(data){
  RDS.data <- as.rds.data.frame(data, id = "RDS_CODE"
                                , recruiter.id = "rec.id"
                                , network.size = "degree")
  
  RDS.data$n.rec <- get.number.of.recruits(RDS.data)
  RDS.data$RDS.degree <- ifelse(RDS.data$rec.id !="0", RDS.data$n.rec+1, RDS.data$n.rec)
  RDS.data$degree1 <- apply(cbind(RDS.data$RDS.degree, RDS.data$degree), 1, max)
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

RDS.IT <- get.RDS.data(sialon.IT)
RDS.LT <- get.RDS.data(sialon.LT)
RDS.RO <- get.RDS.data(sialon.RO)
RDS.SK <- get.RDS.data(sialon.SK)

print("2.done")


sialon.est <- function(data, RDS.data){
  res <- data.frame()
  aa <- data %>% 
    summarise(n = n(), sum(HSU, na.rm = TRUE), mean(age, na.rm = TRUE),  mean(HIV, na.rm = TRUE), mean(HSU, na.rm = TRUE)
              , mean(HIV_testing_history, na.rm = TRUE)
              , sum(ART_coverage[which(HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
              , mean(homosexual, na.rm = TRUE)
              , mean(other_testing_history, na.rm = TRUE)
              , mean(howmany_nonsteady_malepartners, na.rm = TRUE)
              , mean(howmany_nonsteady_malepartners_unprotected, na.rm = TRUE)
              , mean(injected_drug, na.rm = TRUE))

  yourself <- c()
  for(j in 1:7){
    yourself[j] <- length(which(data$yourself==j))/(dim(data)[1]-sum(is.na(data$yourself)))
  }
  occupation <- c()
  for(j in 1:7){
    occupation[j] <- length(which(data$current_occupation==j))/(dim(data)[1]-sum(is.na(data$current_occupation)))
  }
  
  res <- as.data.frame(cbind(aa, t(yourself), t(occupation)))
  names(res) <- c("n", "n.HSU", "age", "HIV", "HSU", "HIV testing history"
                  , "ART", "homosexual", "othertesting history"
                  , "howmany_nonsteady_malepartners", "howmany_nonsteady_malepartners_unprotected"
                  , "injected_drug", "y1.homosexual", "y2.bisexual", "y3.stright"
                  , "y4.any other", "y5.don't use a term", "y6.don't know", "y7.other"
                  , "o1.employed", "o2.self employed", "o3.unemployed", "o4.student"
                  , "o5.retired", "o6.sick leave", "o7.other")
  
  myvar <- c("age", "HIV", "HSU", "HIV_testing_history"
             , "ART_coverage", "homosexual", "other_testing_history"
             , "howmany_nonsteady_malepartners", "howmany_nonsteady_malepartners_unprotected"
             , "injected_drug")
  for(k in 1:length(myvar)){
    if(myvar[k] == "ART_coverage"){
      assign(myvar[k], RDS.II.estimates(rds.data=RDS.data ,outcome.variable=myvar[k], subset = HIV == 1))
    }else{
      assign(myvar[k], RDS.II.estimates(rds.data=RDS.data ,outcome.variable=myvar[k]))
    }
  }
  
  yourself.var <- c("y1.homosexual", "y2.bisexual", "y3.stright"
                    , "y4.any other", "y5.don't use a term", "y6.don't know", "y7.other")
  for(j in 1:7){
    RDS.data$yourself.II <- ifelse(RDS.data$yourself==j, 1, 0)
    RDS.data$yourself.II[which(is.na(RDS.data$yourself))] <- NA
    assign(yourself.var[j], RDS.II.estimates(rds.data=RDS.data ,outcome.variable = "yourself.II"))
  }
  
  occupation.var <- c("o1.employed", "o2.self employed", "o3.unemployed", "o4.student"
                    , "o5.retired", "o6.sick leave", "o7.other") 
  
  for(j in 1:7){
    RDS.data$occupation.II <- ifelse(RDS.data$current_occupation==j, 1, 0)
    RDS.data$occupation.II[which(is.na(RDS.data$current_occupation))] <- NA
    assign(occupation.var[j], RDS.II.estimates(rds.data=RDS.data ,outcome.variable = "occupation.II"))
  }
  
  
  var.RDS <- c(myvar,yourself.var, occupation.var)
  res.RDS <- c()
  for(k in 1:length(var.RDS)){
   res.RDS <-  as.data.frame(cbind(res.RDS, get(var.RDS[k])$interval[1,1:3]))
  }
  
  colnames(res.RDS) <- var.RDS
  #colnames(res.RDS) <- c("Estimate", "95% LCI", "95% UCI")
  
  
  return(list(res = res, res.RDS = res.RDS))
}


est.IT <- sialon.est(sialon.IT, RDS.IT)
est.LT <- sialon.est(sialon.LT, RDS.LT)
est.RO <- sialon.est(sialon.RO, RDS.RO)
est.SK <- sialon.est(sialon.SK, RDS.SK)


print("3.done")





########################################################################
## Snowball
########################################################################
snowball.data <- function(RDS.data, seed.var){
  #forsample.data <- data.frame(RDS.data[,c(1,6,22,23,24,11,19,20,21)])
  forsample.data <- data.frame(RDS.data[,c(2,7,25,26,27,12,22,23,24)])
  snow.data <- data.frame()
  while(dim(snow.data)[1] < dim(RDS.data)[1]){
    seed.tmp <- forsample.data[sample(which(forsample.data[[seed.var]]==1),1),]
    snow.data <- rbind(snow.data,data.frame(seed = "0", snow.sample = seed.tmp$RDS_CODE))
    #print(seed.tmp$degree1)
    if(seed.tmp$RDS.degree==1 & seed.tmp$rec.id != "0"){
      snow.sample <- sample(seed.tmp$rec.id, seed.tmp$degree1, replace = TRUE)
    }else if(seed.tmp$RDS.degree==2 & seed.tmp$rec.id != "0"){
      snow.sample <- sample(c(seed.tmp$rec.id,seed.tmp$C1), seed.tmp$degree1, replace = TRUE)
    }else if(seed.tmp$RDS.degree==3 & seed.tmp$rec.id != "0"){
      snow.sample <- sample(c(seed.tmp$rec.id,seed.tmp$C1, seed.tmp$C2), seed.tmp$degree1, replace = TRUE)
    }else if(seed.tmp$RDS.degree==4 & seed.tmp$rec.id != "0"){
      snow.sample <- sample(c(seed.tmp$rec.id,seed.tmp$C1, seed.tmp$C2,seed.tmp$C3), seed.tmp$degree1, replace = TRUE)
    }else if(seed.tmp$RDS.degree==1 & seed.tmp$rec.id == "0"){
      snow.sample <- sample(seed.tmp$C1, seed.tmp$degree1, replace = TRUE)
    }else if(seed.tmp$RDS.degree==2 & seed.tmp$rec.id == "0"){
      snow.sample <- sample(c(seed.tmp$C1, seed.tmp$C2), seed.tmp$degree1, replace = TRUE)
    }else if(seed.tmp$RDS.degree==3 & seed.tmp$rec.id == "0"){
      snow.sample <- sample(c(seed.tmp$C1, seed.tmp$C2, seed.tmp$C3), seed.tmp$degree1, replace = TRUE)
    }
    snow.data <- rbind(snow.data,data.frame(seed = rep(seed.tmp$RDS_CODE, seed.tmp$degree1), snow.sample))
  }
  snow.data <- snow.data[1:dim(RDS.data)[1],]
  
  snow.data$RDS_CODE <- snow.data$snow.sample
  snow.withcov <- left_join(snow.data, RDS.data, by = "RDS_CODE")
  
  snow.withcov$snow_ID <- as.character(1:dim(RDS.data)[1])
  snow.withcov$snow_recid <- NA
  snow.withcov$snow_recid[which(snow.withcov$seed == 0)] <- "0"
  snow.recid <- which(snow.withcov$seed == 0)
  if(dim(snow.withcov)[1] %in% snow.recid){
    snow.recid <- snow.recid[1:(length(snow.recid)-1)]
  }
  for(i in 1:length(snow.recid)){
    if(i != length(snow.recid)){
      snow.withcov$snow_recid[(snow.recid[i]+1):(snow.recid[(i+1)]-1)] <- as.character(snow.recid[i])
    }else{
      snow.withcov$snow_recid[(snow.recid[i]+1):dim(snow.withcov)[1]] <- as.character(snow.recid[i])
    }
  }
  return(snow.withcov)
}


snowball.all.HSU <- function(RDS.data, seed.var){
  seed <- which(RDS.data[[seed.var]]==1)
  snow.data <- data.frame()
  snow.data <- rbind(snow.data,data.frame(seed = rep("0", length(seed)), snow.sample = RDS.data$RDS_CODE[seed]))
  
  tmp.sample <- cbind(RDS.data$rec.id, RDS.data$C1, RDS.data$C2, RDS.data$C3)
  tmp.sample[which(tmp.sample[,1]=="0"),1] <- NA
  rec.id <- seed
  
  while(dim(snow.data)[1]<dim(RDS.data)[1]){
    snow.sample <- rep(NA, length(rec.id))
    for(i in 1:length(rec.id)){
      snow.sample[i] <- sample(na.omit(tmp.sample[rec.id[i],]),1)
    }
    snow.data <- rbind(snow.data,data.frame(seed = RDS.data$RDS_CODE[rec.id], snow.sample))
    #rec.id <- match(snow.sample, RDS.data$RDS_CODE)
  }
  
  snow.data <- snow.data[1:dim(RDS.data)[1],]
  
  snow.data$RDS_CODE <- snow.data$snow.sample
  snow.withcov <- left_join(snow.data, RDS.data, by = "RDS_CODE")
  
  snow.withcov$snow_ID <- as.character(1:dim(RDS.data)[1])
  snow.withcov$snow_recid <- NA
  snow.withcov$snow_recid[which(snow.withcov$seed == 0)] <- "0"
  snow.withcov$snow_recid[(length(seed)+1):dim(RDS.data)[1]] <- as.character(rep(1:length(seed), ceiling(dim(RDS.data)[1]/length(seed)))[(length(seed)+1):dim(RDS.data)[1]])
  return(snow.withcov)
}

snowball.wave1 <- function(RDS.data, sample.size = dim(RDS.data)[1]
                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var){
  
  tmp.sample <- cbind(RDS.data$RDS_CODE, RDS.data$rec.id, RDS.data$C1, RDS.data$C2, RDS.data$C3)
  tmp.sample[which(tmp.sample[,2]=="0"),2] <- NA
  colnames(tmp.sample) <- c("RDS_CODE", "rec.id", "C1", "C2", "C3")
  
  snow.data <- data.frame()
  i <- 1
  k <- 1
  while(dim(snow.data)[1]<dim(RDS.data)[1]){
    seed.tmp <- sample(which(RDS.data[[seed.var]]==1),1)
    snow.data <- rbind(snow.data, data.frame(RDS.recid = "0", RDS.sample.ID = RDS.data$RDS_CODE[seed.tmp]
                                             , wave = 0, RDS.new.recid = "0", RDS.new.ID = as.character(k)
                                             , rec.num = sample(0:number.of.coupons, 1)))
    
    if(snow.data$rec.num[i] == 0){
      k <- k+1
    }else{
      aa <- which(snow.data$RDS.sample.ID[i] == tmp.sample[,1])
      rec.sample <- sample(na.omit(tmp.sample[aa,2:5]), snow.data$rec.num[i], replace = TRUE)
      snow.data <- rbind(snow.data, data.frame(RDS.recid = rep(tmp.sample[aa,1], snow.data$rec.num[i])
                                               , RDS.sample.ID = rec.sample, wave = 1
                                               , RDS.new.recid = rep(snow.data$RDS.new.ID[i], snow.data$rec.num[i])
                                               , RDS.new.ID = as.character((k+1):(k+snow.data$rec.num[i]))
                                               , rec.num = sample(0:number.of.coupons, snow.data$rec.num[i], replace = TRUE)))
      k <- k+snow.data$rec.num[i]+1
    }
    i <- k
  }
  
  snow.data <- snow.data[1:dim(RDS.data)[1],]
  rownames(snow.data) <- as.character(1:dim(RDS.data)[1])
  
  snow.data$RDS_CODE <- snow.data$RDS.sample.ID
  snow.withcov <- left_join(snow.data, RDS.data, by = "RDS_CODE")
  
  #snow.withcov$RDS_ID <- as.character(1:dim(RDS.data)[1])
  #snow.withcov$RDS_recid <- NA
  #snow.withcov$RDS_recid[which(snow.withcov$seed == 0)] <- "0"
  #snow.withcov$snow_recid[(length(seed)+1):dim(RDS.data)[1]] <- as.character(rep(1:length(seed), ceiling(dim(RDS.data)[1]/length(seed)))[(length(seed)+1):dim(RDS.data)[1]])
  return(snow.withcov)
}





snowball.wave2 <- function(RDS.data, sample.size = dim(RDS.data)[1]
                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var){
  
  tmp.sample <- cbind(RDS.data$RDS_CODE, RDS.data$rec.id, RDS.data$C1, RDS.data$C2, RDS.data$C3)
  tmp.sample[which(tmp.sample[,2]=="0"),2] <- NA
  colnames(tmp.sample) <- c("RDS_CODE", "rec.id", "C1", "C2", "C3")
  
  snow.data <- data.frame()
  i <- 1
  k <- 1
  while(dim(snow.data)[1]<dim(RDS.data)[1]){
    seed.tmp <- sample(which(RDS.data[[seed.var]]==1),1)
    snow.data <- rbind(snow.data, data.frame(RDS.recid = "0", RDS.sample.ID = RDS.data$RDS_CODE[seed.tmp]
                                             , wave = 0, RDS.new.recid = "0", RDS.new.ID = as.character(k)
                                             , rec.num = sample(0:number.of.coupons, 1)))
    
    if(snow.data$rec.num[i] == 0){
      k <- k+1
    }else{
      aa <- which(snow.data$RDS.sample.ID[i] == tmp.sample[,1])
      rec.sample <- sample(na.omit(tmp.sample[aa,2:5]), snow.data$rec.num[i], replace = TRUE)
      snow.data <- rbind(snow.data, data.frame(RDS.recid = rep(tmp.sample[aa,1], snow.data$rec.num[i])
                                               , RDS.sample.ID = rec.sample, wave = 1
                                               , RDS.new.recid = rep(snow.data$RDS.new.ID[i], snow.data$rec.num[i])
                                               , RDS.new.ID = as.character((k+1):(k+snow.data$rec.num[i]))
                                               , rec.num = sample(0:number.of.coupons, snow.data$rec.num[i], replace = TRUE)))
      k <- k+snow.data$rec.num[i]
      
      
      for(j in (i+1):(i+snow.data$rec.num[i]))
        if(snow.data$rec.num[j] == 0){
          # k <- k+1
        }else{
          aa <- which(snow.data$RDS.sample.ID[j] == tmp.sample[,1])
          rec.sample <- sample(na.omit(tmp.sample[aa,2:5]), snow.data$rec.num[j], replace = TRUE)
          snow.data <- rbind(snow.data, data.frame(RDS.recid = rep(tmp.sample[aa,1], snow.data$rec.num[j])
                                                   , RDS.sample.ID = rec.sample, wave = 2
                                                   , RDS.new.recid = rep(snow.data$RDS.new.ID[j], snow.data$rec.num[j])
                                                   , RDS.new.ID = as.character((k+1):(k+snow.data$rec.num[j]))
                                                   , rec.num = 0))
          
          k <- k+snow.data$rec.num[j]
        }
      k <- k+1
    }
    i <- k
  }
  
  snow.data <- snow.data[1:dim(RDS.data)[1],]
  rownames(snow.data) <- as.character(1:dim(RDS.data)[1])
  
  snow.data$RDS_CODE <- snow.data$RDS.sample.ID
  snow.withcov <- left_join(snow.data, RDS.data, by = "RDS_CODE")
  
  #snow.withcov$RDS_ID <- as.character(1:dim(RDS.data)[1])
  #snow.withcov$RDS_recid <- NA
  #snow.withcov$RDS_recid[which(snow.withcov$seed == 0)] <- "0"
  #snow.withcov$snow_recid[(length(seed)+1):dim(RDS.data)[1]] <- as.character(rep(1:length(seed), ceiling(dim(RDS.data)[1]/length(seed)))[(length(seed)+1):dim(RDS.data)[1]])
  return(snow.withcov)
}

snowball.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway, seed.var = seed.var){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data, seed.var = seed.var)
  
  aa <- snow.data %>% 
    summarise(n = n(), mean(age, na.rm = TRUE), mean(HSU, na.rm = TRUE), mean(HIV, na.rm = TRUE)
              , sum(HIV, na.rm = TRUE), mean(HIV_testing_history, na.rm = TRUE)
              , sum(ART_coverage[which(snow.data$HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
              , mean(homosexual, na.rm = TRUE))
  bb <- c()
  for(j in 1:7){
    bb[j] <- length(which(snow.data$yourself==j))/dim(snow.data)[1]
  }
  HSU <- data.frame()
  
  snow.data1 <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
  
  a00 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HSU')
  a0 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV')
  a1 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV_testing_history')
  if(length(which(snow.data1$HIV==1))==0){
    a2 <- a1
    a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
  }else{
    a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
  }
  a3 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='homosexual')
  a4 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='age')
  HSU <- as.data.frame(cbind(aa, a00$estimate, a00$interval[2], a00$interval[3]
                             , a0$estimate, a0$interval[2], a0$interval[3]
                             , a1$estimate, a1$interval[2], a1$interval[3]
                             , a2$estimate, a2$interval[2], a2$interval[3]
                             , t(bb), a3$estimate, a3$interval[2], a3$interval[3]
                             , a4$estimate, a4$interval[2], a4$interval[3]))
  
  
  
  for(i in 2:iter){
    set.seed(setseed+i)
    snow.data <- snowball.sampleway(RDS.data, seed.var = seed.var)
    
    aa <- snow.data %>% 
      summarise(n = n(), mean(age, na.rm = TRUE), mean(HSU, na.rm = TRUE),  mean(HIV, na.rm = TRUE)
                , sum(HIV, na.rm = TRUE), mean(HIV_testing_history, na.rm = TRUE)
                , sum(ART_coverage[which(snow.data$HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
                , mean(homosexual, na.rm = TRUE))
    bb <- c()
    for(j in 1:7){
      bb[j] <- length(which(snow.data$yourself==j))/dim(snow.data)[1]
    }
    
    snow.data1 <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                    , recruiter.id = "RDS.new.recid"
                                    , network.size = "degree")
    
    a00 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HSU')
    a0 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV')
    a1 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV_testing_history')
    if(length(which(snow.data1$HIV==1))==0){
      a2 <- a1
      a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
    }else{
      a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
    }
    a3 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='homosexual')
    a4 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='age')
    HSU <- rbind(HSU, as.data.frame(cbind(aa, a00$estimate, a00$interval[2], a00$interval[3]
                                          , a0$estimate, a0$interval[2], a0$interval[3]
                                          , a1$estimate, a1$interval[2], a1$interval[3]
                                          , a2$estimate, a2$interval[2], a2$interval[3]
                                          , t(bb), a3$estimate, a3$interval[2], a3$interval[3]
                                          , a4$estimate, a4$interval[2], a4$interval[3])))
    
    
    print(i)
  }
  return(HSU) 
}


if(2==3){


snowball.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data)
  
  aa <- snow.data %>% 
    summarise(n = n(), mean(age, na.rm = TRUE), mean(HSU, na.rm = TRUE), mean(HIV, na.rm = TRUE)
              , sum(HIV, na.rm = TRUE), mean(HIV_testing_history, na.rm = TRUE)
              , sum(ART_coverage[which(snow.data$HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
              , mean(homosexual, na.rm = TRUE))
  bb <- c()
  for(j in 1:7){
    bb[j] <- length(which(snow.data$yourself==j))/dim(snow.data)[1]
  }
  HSU <- data.frame()
  
  snow.data1 <- as.rds.data.frame(snow.data, id = "snow_ID"
                                  , recruiter.id = "snow_recid"
                                  , network.size = "degree")
  
  a00 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HSU')
  a0 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV')
  a1 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV_testing_history')
  if(length(which(snow.data1$HIV==1))==0){
    a2 <- a1
    a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
  }else{
    a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
  }
  a3 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='homosexual')
  a4 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='age')
  HSU <- as.data.frame(cbind(aa, a00$estimate, a00$interval[2], a00$interval[3]
                             , a0$estimate, a0$interval[2], a0$interval[3]
                             , a1$estimate, a1$interval[2], a1$interval[3]
                             , a2$estimate, a2$interval[2], a2$interval[3]
                             , t(bb), a3$estimate, a3$interval[2], a3$interval[3]
                             , a4$estimate, a4$interval[2], a4$interval[3]))
  
  
  
  for(i in 2:iter){
    set.seed(setseed+i)
    snow.data <- snowball.sampleway(RDS.data)
    
    aa <- snow.data %>% 
      summarise(n = n(), mean(age, na.rm = TRUE), mean(HSU, na.rm = TRUE),  mean(HIV, na.rm = TRUE)
                , sum(HIV, na.rm = TRUE), mean(HIV_testing_history, na.rm = TRUE)
                , sum(ART_coverage[which(snow.data$HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
                , mean(homosexual, na.rm = TRUE))
    bb <- c()
    for(j in 1:7){
      bb[j] <- length(which(snow.data$yourself==j))/dim(snow.data)[1]
    }
    
    snow.data1 <- as.rds.data.frame(snow.data, id = "snow_ID"
                                    , recruiter.id = "snow_recid"
                                    , network.size = "degree")
    
    a00 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HSU')
    a0 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV')
    a1 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV_testing_history')
    if(length(which(snow.data1$HIV==1))==0){
      a2 <- a1
      a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
    }else{
      a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
    }
    a3 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='homosexual')
    a4 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='age')
    HSU <- rbind(HSU, as.data.frame(cbind(aa, a00$estimate, a00$interval[2], a00$interval[3]
                                          , a0$estimate, a0$interval[2], a0$interval[3]
                                          , a1$estimate, a1$interval[2], a1$interval[3]
                                          , a2$estimate, a2$interval[2], a2$interval[3]
                                          , t(bb), a3$estimate, a3$interval[2], a3$interval[3]
                                          , a4$estimate, a4$interval[2], a4$interval[3])))
    
    
    print(i)
  }
  return(HSU) 
}


IT.HSU.rep <- snowball.rep(RDS.data = RDS.IT, 1001, iter = 1000, snowball.sampleway = snowball.data)
LT.HSU.rep <- snowball.rep(RDS.data = RDS.LT, 1001, iter = 1000, snowball.sampleway = snowball.data)
RO.HSU.rep <- snowball.rep(RDS.data = RDS.RO, 1001, iter = 1000, snowball.sampleway = snowball.data)
SK.HSU.rep <- snowball.rep(RDS.data = RDS.SK, 1001, iter = 1000, snowball.sampleway = snowball.data)

print("4.done")

save.image("sialon_all_rep.RData")



IT.all.rep <- snowball.rep(RDS.data = RDS.IT, 1001, iter = 1000, snowball.sampleway = snowball.all.HSU)
LT.all.rep <- snowball.rep(RDS.data = RDS.LT, 1001, iter = 1000, snowball.sampleway = snowball.all.HSU)
RO.all.rep <- snowball.rep(RDS.data = RDS.RO, 1001, iter = 1000, snowball.sampleway = snowball.all.HSU)
SK.all.rep <- snowball.rep(RDS.data = RDS.SK, 1001, iter = 1000, snowball.sampleway = snowball.all.HSU)


print("5.done")
save.image("sialon_all_rep.RData")
}


IT.HSU.wave1.rep <- snowball.rep(RDS.data = RDS.IT, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HSU")
LT.HSU.wave1.rep <- snowball.rep(RDS.data = RDS.LT, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HSU")
RO.HSU.wave1.rep <- snowball.rep(RDS.data = RDS.RO, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HSU")
SK.HSU.wave1.rep <- snowball.rep(RDS.data = RDS.SK, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HSU")
save.image("sialon_all_rep.RData")

IT.HIV.wave1.rep <- snowball.rep(RDS.data = RDS.IT, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HIV_testing_history")
LT.HIV.wave1.rep <- snowball.rep(RDS.data = RDS.LT, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HIV_testing_history")
RO.HIV.wave1.rep <- snowball.rep(RDS.data = RDS.RO, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HIV_testing_history")
SK.HIV.wave1.rep <- snowball.rep(RDS.data = RDS.SK, 1001, iter = 1000, snowball.sampleway = snowball.wave1, seed.var = "HIV_testing_history")
save.image("sialon_all_rep.RData")

IT.HSU.wave2.rep <- snowball.rep(RDS.data = RDS.IT, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HSU")
LT.HSU.wave2.rep <- snowball.rep(RDS.data = RDS.LT, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HSU")
RO.HSU.wave2.rep <- snowball.rep(RDS.data = RDS.RO, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HSU")
SK.HSU.wave2.rep <- snowball.rep(RDS.data = RDS.SK, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HSU")
save.image("sialon_all_rep.RData")

IT.HIV.wave2.rep <- snowball.rep(RDS.data = RDS.IT, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HIV_testing_history")
LT.HIV.wave2.rep <- snowball.rep(RDS.data = RDS.LT, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HIV_testing_history")
RO.HIV.wave2.rep <- snowball.rep(RDS.data = RDS.RO, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HIV_testing_history")
SK.HIV.wave2.rep <- snowball.rep(RDS.data = RDS.SK, 1001, iter = 1000, snowball.sampleway = snowball.wave2, seed.var = "HIV_testing_history")
save.image("sialon_all_rep.RData")


########################################################################
## RDS sample
########################################################################
RDS.sample.dongah <- function(RDS.data, sample.size = dim(RDS.data)[1]
                              , number.of.coupons = 3, sample.with.replacement = TRUE){
  number.of.seeds <- length(unique(get.seed.id(RDS.data)))
  tmp.sample <- cbind(RDS.data$RDS_CODE, RDS.data$rec.id, RDS.data$C1, RDS.data$C2, RDS.data$C3)
  tmp.sample[which(tmp.sample[,2]=="0"),2] <- NA
  colnames(tmp.sample) <- c("RDS_CODE", "rec.id", "C1", "C2", "C3")
  
  
  seed <- sample(1:dim(RDS.data)[1], number.of.seeds)
  snow.data <- data.frame(RDS.recid = "0", RDS.sample.ID = RDS.data$RDS_CODE[seed], wave = 0, RDS.new.recid = "0", RDS.new.ID = as.character(1:number.of.seeds))
  j <- 1  ## wave 
  k <- 1  ## n.participants in j-1 waves
  l <- number.of.seeds+1  ## new.ID
  
  #rec.num1 <- sample(0:number.of.coupons, dim(RDS.data)[1], replace = TRUE)
  while(dim(snow.data)[1]<dim(RDS.data)[1]){
    for(i in k:(k+length(which(snow.data$wave==(j-1)))-1)){
      rec.num <- sample(0:number.of.coupons, 1, prob = c(1.3,2,2,2))
      #print(rec.num)
      if(rec.num==0){
        q <- 1
      }else{
        aa <- which(snow.data$RDS.sample.ID[i] == tmp.sample[,1])
        rec.sample <- sample(na.omit(tmp.sample[aa,2:5]), rec.num, replace = TRUE)
        snow.data <- rbind(snow.data, data.frame(RDS.recid = rep(tmp.sample[aa,1], rec.num)
                                                 , RDS.sample.ID = rec.sample, wave = rep(j, rec.num)
                                                 , RDS.new.recid = rep(snow.data$RDS.new.ID[i], rec.num)
                                                 , RDS.new.ID = as.character(l:(l+rec.num-1))))
        l <- l+rec.num ## RDS.new.id
        
      }
      
    }
    k <- k+length(which(snow.data$wave==(j-1))) ## n.participants in j-1 waves
    j <- j+1 ## wave
  }
  snow.data <- snow.data[1:dim(RDS.data)[1],]
  rownames(snow.data) <- as.character(1:dim(RDS.data)[1])
  
  snow.data$RDS_CODE <- snow.data$RDS.sample.ID
  snow.withcov <- left_join(snow.data, RDS.data, by = "RDS_CODE")
  
  #snow.withcov$RDS_ID <- as.character(1:dim(RDS.data)[1])
  #snow.withcov$RDS_recid <- NA
  #snow.withcov$RDS_recid[which(snow.withcov$seed == 0)] <- "0"
  #snow.withcov$snow_recid[(length(seed)+1):dim(RDS.data)[1]] <- as.character(rep(1:length(seed), ceiling(dim(RDS.data)[1]/length(seed)))[(length(seed)+1):dim(RDS.data)[1]])
  return(snow.withcov)
  
  
}

snowball.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data)
  
  aa <- snow.data %>% 
    summarise(n = n(), mean(age, na.rm = TRUE), mean(HSU, na.rm = TRUE),  mean(HIV, na.rm = TRUE)
              , sum(HIV, na.rm = TRUE), mean(HIV_testing_history, na.rm = TRUE)
              , sum(ART_coverage[which(snow.data$HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
              , mean(homosexual, na.rm = TRUE))
  bb <- c()
  for(j in 1:7){
    bb[j] <- length(which(snow.data$yourself==j))/dim(snow.data)[1]
  }
  
  HSU <- data.frame()
  
  snow.data1 <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
  
  a00 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HSU')
  a0 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV')
  a1 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV_testing_history')
  if(length(which(snow.data1$HIV==1))==0){
    a2 <- a1
    a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
  }else{
    a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
  }
  a3 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='homosexual')
  a4 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='age')
  
  #SPRTBA.fix_test.per_hsu <- SPRTBA(net = snow.data1, fix.variable = "HIV_testing_history", permute.variable = "HSU", test = chisq.test, n = 1000, correct = TRUE)
  #SPRTBA.fix_hsu.per_test <- SPRTBA(net = snow.data1, fix.variable = "HSU", permute.variable = "HIV_testing_history", test = chisq.test, n = 1000, correct = TRUE)
  #SPRTBA.fix_hiv.per_hsu <- SPRTBA(net = snow.data1, fix.variable = "HIV", permute.variable = "HSU", test = chisq.test, n = 1000, correct = TRUE)
  #SPRTBA.fix_hsu.per_hiv <- SPRTBA(net = snow.data1, fix.variable = "HSU", permute.variable = "HIV", test = chisq.test, n = 1000, correct = TRUE)
  
  HSU <- as.data.frame(cbind(aa, a00$estimate, a00$interval[2], a00$interval[3]
                             , a0$estimate, a0$interval[2], a0$interval[3]
                             , a1$estimate, a1$interval[2], a1$interval[3]
                             , a2$estimate, a2$interval[2], a2$interval[3]
                             , t(bb), a3$estimate, a3$interval[2], a3$interval[3]
                             , a4$estimate, a4$interval[2], a4$interval[3]))
                            # , SPRTBA.fix_test.per_hsu$p.value.mc, SPRTBA.fix_hsu.per_test$p.value.mc
                            # , SPRTBA.fix_hiv.per_hsu$p.value.mc, SPRTBA.fix_hsu.per_hiv$p.value.mc))
  
  
  for(p in 2:iter){
    if(p==220){
      p <- 1000
    }
    set.seed(setseed+p)
    
    snow.data <- snowball.sampleway(RDS.data)
    
    aa <- snow.data %>% 
      summarise(n = n(), mean(age, na.rm = TRUE), mean(HSU, na.rm = TRUE),  mean(HIV, na.rm = TRUE)
                , sum(HIV, na.rm = TRUE), mean(HIV_testing_history, na.rm = TRUE)
                , sum(ART_coverage[which(snow.data$HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
                , mean(homosexual, na.rm = TRUE))
    bb <- c()
    for(j in 1:7){
      bb[j] <- length(which(snow.data$yourself==j))/dim(snow.data)[1]
    }
    
    snow.data1 <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                    , recruiter.id = "RDS.new.recid"
                                    , network.size = "degree")
    
    a00 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HSU')
    a0 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV')
    a1 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='HIV_testing_history')
    if(length(which(snow.data1$HIV==1))==0){
      a2 <- a1
      a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
    }else{
      a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
    }
    a3 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='homosexual')
    a4 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='age')
    #SPRTBA.fix_test.per_hsu <- SPRTBA(net = snow.data1, fix.variable = "HIV_testing_history", permute.variable = "HSU", test = chisq.test, n = 1000, correct = TRUE)
    #SPRTBA.fix_hsu.per_test <- SPRTBA(net = snow.data1, fix.variable = "HSU", permute.variable = "HIV_testing_history", test = chisq.test, n = 1000, correct = TRUE)
    #SPRTBA.fix_hiv.per_hsu <- SPRTBA(net = snow.data1, fix.variable = "HIV", permute.variable = "HSU", test = chisq.test, n = 1000, correct = TRUE)
    #SPRTBA.fix_hsu.per_hiv <- SPRTBA(net = snow.data1, fix.variable = "HSU", permute.variable = "HIV", test = chisq.test, n = 1000, correct = TRUE)
    
    HSU <- rbind(HSU, as.data.frame(cbind(aa, a00$estimate, a00$interval[2], a00$interval[3]
                                          , a0$estimate, a0$interval[2], a0$interval[3]
                                          , a1$estimate, a1$interval[2], a1$interval[3]
                                          , a2$estimate, a2$interval[2], a2$interval[3]
                                          , t(bb), a3$estimate, a3$interval[2], a3$interval[3]
                                          , a4$estimate, a4$interval[2], a4$interval[3])))
                                          #, SPRTBA.fix_test.per_hsu$p.value.mc, SPRTBA.fix_hsu.per_test$p.value.mc
                                          #, SPRTBA.fix_hiv.per_hsu$p.value.mc, SPRTBA.fix_hsu.per_hiv$p.value.mc)))
    
    print(p)
  }
  return(HSU) 
}



IT.RDS.rep <- snowball.rep(RDS.data = RDS.IT, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)
LT.RDS.rep <- snowball.rep(RDS.data = RDS.LT, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)
RO.RDS.rep <- snowball.rep(RDS.data = RDS.RO, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)
SK.RDS.rep <- snowball.rep(RDS.data = RDS.SK, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)

print("6.done")

save.image("sialon_all_rep.RData")
save.image("RDS_rep.RData")





if(2==3){


head(IT.rep)



#pdf("snowball.HIV.pdf")

IT.rep <- IT.all.rep 
LT.rep <- LT.all.rep 
RO.rep <- RO.all.rep 
SK.rep <- SK.all.rep 
  
makeboxplot <- function(rep.data, mymain, n.est , reach = NULL){
  boxplot(rep.data, col = mycol, xaxt = "n", main = mymain, cex.main = 2, cex.axis = 1.5, ylim = c(0,1))
  axis(1, at=1:5,labels=c("HSU", "HIV", "HIV test history", "ART", "Homosexual"), cex.axis = 1.2)
  segments(x0 = x0s, x1 = x1s, y0 = c(est.res[n.est,4], est.res[n.est,5], est.res[n.est,6], est.res[n.est,7], est.res[n.est,8]), col = "red", lwd = 2)
  if(is.null(reach)==FALSE){
    segments(x0 = x0s, x1 = x1s, y0 = as.numeric(reach[4:8]), col = "blue", lwd = 2, lty = 2)
    segments(x0 = x0s, x1 = x1s, y0 = as.numeric(reach[12:16]), col = "green", lwd = 2, lty = 2)
    legend("topleft", lwd = 2, lty = c(1,2,2), col = c("red", "blue", "green"), legend = c("original RDS data", "reachable", "not reachable"), cex = 1.5)
  }
  legend("topleft", lwd = 2, lty = 1, col = "red", legend = "original RDS data", cex = 1.5)
}

n <- 4
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
mycol <- brewer.pal(n = 5, name = "RdBu")



age <- data.frame(Italy=IT.rep[,2], Lithuania=LT.rep[,2], Romania=RO.rep[,2], Slovakia=SK.rep[,2])
age.RDSII <- data.frame(Italy=IT.rep[,31], Lithuania=LT.rep[,31], Romania=RO.rep[,31], Slovakia=SK.rep[,31])

age.original <- est.res[,c(3,31)]
n <- 5
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
mycol <- brewer.pal(n = 5, name = "RdBu")


pdf("snowball.reachable.rep.prop.pdf")
makeboxplot(IT.rep[,c(3,4,6,7,8)], mymain = "Italy (mean)", n.est = 1, reach = ITreach)
makeboxplot(LT.rep[,c(3,4,6,7,8)], mymain = "Lithuania (mean)", n.est = 2, reach = LTreach)
makeboxplot(RO.rep[,c(3,4,6,7,8)], mymain = "Romania (mean)", n.est = 3, reach = ROreach)
makeboxplot(SK.rep[,c(3,4,6,7,8)], mymain = "Slovakia (mean)", n.est = 4, reach = SKreach)
makeboxplot(IT.rep[,c(9,12,15,18,28)], mymain = "Italy (RDS II est)", n.est = 1)
makeboxplot(LT.rep[,c(9,12,15,18,28)], mymain = "Lithuania (RDS II est)", n.est = 2)
makeboxplot(RO.rep[,c(9,12,15,18,28)], mymain = "Romania (RDS II est)", n.est = 3)
makeboxplot(SK.rep[,c(9,12,15,18,28)], mymain = "Slovakia (RDS II est)", n.est = 4)

n <- 4
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4


boxplot(age, col = mycol, main = "age (mean)", cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
segments(x0 = x0s, x1 = x1s, y0 = age.original[,1], col = "red", lwd = 2)
legend("topleft", lwd = 2, lty = 1, col = "red", legend = "original RDS data", cex = 1.5)

boxplot(age.RDSII, col = mycol, main = "age (RDS II est)", cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
segments(x0 = x0s, x1 = x1s, y0 = age.original[,2], col = "red", lwd = 2)
legend("topleft", lwd = 2, lty = 1, col = "red", legend = "original RDS data", cex = 1.5)

dev.off()



### possible sample wave 1

reach.notreach <- function(reach, notreach){
  res <- data.frame()
  aa <- reach %>% 
    summarise(n = n(), n()/(dim(possible.RDS)[1]+dim(notpossible.RDS)[1]), mean(age, na.rm = TRUE),  mean(HIV, na.rm = TRUE), mean(HSU, na.rm = TRUE)
              , mean(HIV_testing_history, na.rm = TRUE)
              , sum(ART_coverage[which(HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
              , mean(homosexual, na.rm = TRUE))
  bb <- notreach %>% 
    summarise(n = n(), n()/(dim(possible.RDS)[1]+dim(notpossible.RDS)[1]), mean(age, na.rm = TRUE),  mean(HIV, na.rm = TRUE), mean(HSU, na.rm = TRUE)
              , mean(HIV_testing_history, na.rm = TRUE)
              , sum(ART_coverage[which(HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
              , mean(homosexual, na.rm = TRUE))
  #a00 <- RDS.II.estimates(rds.data=RDS.data ,outcome.variable='HSU')
  #a0 <- RDS.II.estimates(rds.data=RDS.data ,outcome.variable='HIV')
  #a1 <- RDS.II.estimates(rds.data=RDS.data ,outcome.variable='HIV_testing_history')
  #a2 <- RDS.II.estimates(rds.data=RDS.data ,outcome.variable='ART_coverage', subset = HIV == 1)
  
  res <- as.data.frame(cbind(aa, bb))
  #a00$estimate, a00$interval[2], a00$interval[3]
  #                          , a0$estimate, a0$interval[2], a0$interval[3]
  #                         , a1$estimate, a1$interval[2], a1$interval[3]
  #                        , a2$estimate, a2$interval[2], a2$interval[3]))
  
  return(res)
}


RDS.data <- RDS.SK
#possible <- union(union(union(union(RDS.data[which(RDS.data$HSU==1),"RDS_CODE"],RDS.data[which(RDS.data$HSU==1),"C1"]),RDS.data[which(RDS.data$HSU==1),"C2"]),RDS.data[which(RDS.data$HSU==1),"C3"]), RDS.data[which(RDS.data$HSU==1),"rec.id"])
possible <- union(union(union(union(RDS.data[which(RDS.data$HIV_testing_history==1),"RDS_CODE"],RDS.data[which(RDS.data$HIV_testing_history==1),"C1"]),RDS.data[which(RDS.data$HIV_testing_history==1),"C2"]),RDS.data[which(RDS.data$HIV_testing_history==1),"C3"]), RDS.data[which(RDS.data$HIV_testing_history==1),"rec.id"])

notpossible <- setdiff(RDS.data$RDS_CODE, possible)
length(possible)
possible.RDS <- RDS.data[na.omit(match(possible, RDS.data$RDS_CODE)),]
notpossible.RDS <- RDS.data[na.omit(match(notpossible, RDS.data$RDS_CODE)),]
dim(possible.RDS)
dim(notpossible.RDS)

possible.RDS %>% 
  summarise(n = n(), sum(HSU, na.rm = TRUE), mean(age, na.rm = TRUE), mean(HIV, na.rm = TRUE)
            , mean(HIV_testing_history, na.rm = TRUE)
            , sum(ART_coverage[which(HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
            , mean(homosexual, na.rm = TRUE), mean(HSU, na.rm = TRUE), sum(HIV, na.rm = TRUE))
  
notpossible.RDS %>% 
  summarise(n = n(), sum(HSU, na.rm = TRUE), mean(age, na.rm = TRUE), mean(HIV, na.rm = TRUE)
            , mean(HIV_testing_history, na.rm = TRUE)
            , sum(ART_coverage[which(HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
            , mean(homosexual, na.rm = TRUE), mean(HSU, na.rm = TRUE), sum(HIV, na.rm = TRUE))





ITreach <- reach.notreach(possible.RDS, notpossible.RDS)
LTreach <- reach.notreach(possible.RDS, notpossible.RDS)
ROreach <- reach.notreach(possible.RDS, notpossible.RDS)
SKreach <- reach.notreach(possible.RDS, notpossible.RDS)






t.test(possible.RDS$age, notpossible.RDS$age)
t.test(possible.RDS$HIV, notpossible.RDS$HIV)
t.test(possible.RDS$HSU, notpossible.RDS$HSU)

t.test(possible.RDS$HIV_testing_history, notpossible.RDS$HIV_testing_history)
t.test(possible.RDS$ART_coverage, notpossible.RDS$ART_coverage)
t.test(possible.RDS$homosexual, notpossible.RDS$homosexual)



RDS.II.estimates(rds.data=possible.RDS ,outcome.variable='HSU')
RDS.II.estimates(rds.data=possible.RDS ,outcome.variable='HIV')
RDS.II.estimates(rds.data=possible.RDS ,outcome.variable='HIV_testing_history')



reingold.tilford.plot(RDS.RO, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Romania")

RDS.SK %>% 
  summarise(n = n(), sum(HSU, na.rm = TRUE), mean(age, na.rm = TRUE), mean(HIV, na.rm = TRUE)
            , mean(HIV_testing_history, na.rm = TRUE)
            , sum(ART_coverage[which(HIV==1)], na.rm = TRUE)/sum(HIV, na.rm = TRUE)
            , mean(homosexual, na.rm = TRUE), mean(HSU, na.rm = TRUE), sum(HIV, na.rm = TRUE))

RDS.II.estimates(rds.data=RDS.IT ,outcome.variable='age')

HIV_testing <- data.frame(Italy=IT.rep[,6], Lithuania=LT.rep[,6], Romania=RO.rep[,6], Slovakia=SK.rep[,6])
HIV_testing.original <- est.res[,c(6,15)]
HIV_testing.reachable <- as.numeric(c(ITreach[6], LTreach[6], ROreach[6], SKreach[6]))
HIV_testing.notreachable <- as.numeric(c(ITreach[14], LTreach[14], ROreach[14], SKreach[14]))



n <- 4
x0s <- 1:n - 0.4
x1s <- 1:n + 0.4
pdf("HIVhistry.reachable.pdf")
boxplot(HIV_testing, col = mycol, main = "HIV testing history (mean)", cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
segments(x0 = x0s, x1 = x1s, y0 = HIV_testing.original[,1], col = "red", lwd = 2)
segments(x0 = x0s, x1 = x1s, y0 = HIV_testing.reachable, col = "blue", lwd = 2, lty = 2)
segments(x0 = x0s, x1 = x1s, y0 = HIV_testing.notreachable, col = "green", lwd = 2, lty = 2)
legend("topright", lwd = 2, lty = c(1,2,2), col = c("red", "blue", "green"), legend = c("original RDS data", "reachable", "not reachable"), cex = 1.5)
dev.off()

diff.original <- diff.reachable <- array(dim= dim(HIV_testing))
for(i in 1:4){
  diff.original[,i] <- HIV_testing[,i]-HIV_testing.original[i,1]
  diff.reachable[,i] <- HIV_testing[,i]-HIV_testing.reachable[i]
}
plot(diff.reachable[,1], diff.original[,1])
abline(a=0,b=1, col = 2)



### possible sample wave 2
RDS.data <- RDS.SK
#possible1 <- union(union(union(union(RDS.data[which(RDS.data$HSU==1),"RDS_CODE"],RDS.data[which(RDS.data$HSU==1),"C1"]),RDS.data[which(RDS.data$HSU==1),"C2"]),RDS.data[which(RDS.data$HSU==1),"C3"]), RDS.data[which(RDS.data$HSU==1),"rec.id"])
possible1 <- union(union(union(union(RDS.data[which(RDS.data$HIV_testing_history==1),"RDS_CODE"],RDS.data[which(RDS.data$HIV_testing_history==1),"C1"]),RDS.data[which(RDS.data$HIV_testing_history==1),"C2"]),RDS.data[which(RDS.data$HIV_testing_history==1),"C3"]), RDS.data[which(RDS.data$HIV_testing_history==1),"rec.id"])

aa <- union(possible1, na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "C1"]))
bb <- union(aa,na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "C2"]))
cc <- union(bb,na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "C3"]))
possible <- union(cc,na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "rec.id"]))
notpossible <- setdiff(RDS.data$RDS_CODE, possible)
length(possible)
possible.RDS <- RDS.data[na.omit(match(possible, RDS.data$RDS_CODE)),]
notpossible.RDS <- RDS.data[na.omit(match(notpossible, RDS.data$RDS_CODE)),]
dim(possible.RDS)
dim(notpossible.RDS)


ITreach <- reach.notreach(possible.RDS, notpossible.RDS)
LTreach <- reach.notreach(possible.RDS, notpossible.RDS)
ROreach <- reach.notreach(possible.RDS, notpossible.RDS)
SKreach <- reach.notreach(possible.RDS, notpossible.RDS)


t.test(possible.RDS$age, notpossible.RDS$age)
t.test(possible.RDS$HIV, notpossible.RDS$HIV)
t.test(possible.RDS$HSU, notpossible.RDS$HSU)
t.test(possible.RDS$HIV_testing_history, notpossible.RDS$HIV_testing_history)
t.test(possible.RDS$ART_coverage, notpossible.RDS$ART_coverage)
t.test(possible.RDS$homosexual, notpossible.RDS$homosexual)




HSU.wave1 <- snowball.wave1(RDS.data = RDS.SK, sample.size = dim(RDS.data)[1]
               , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HSU")

HSU.wave2 <- snowball.wave2(RDS.data = RDS.SK, sample.size = dim(RDS.data)[1]
                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HSU")
HIV.wave1 <- snowball.wave1(RDS.data = RDS.SK, sample.size = dim(RDS.data)[1]
                            , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HIV_testing_history")

HIV.wave2 <- snowball.wave2(RDS.data = RDS.SK, sample.size = dim(RDS.data)[1]
                            , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HIV_testing_history")

  
reRDS <- RDS.sample.dongah(RDS.data = RDS.SK, sample.size = dim(RDS.data)[1]
                              , number.of.coupons = 3, sample.with.replacement = TRUE)

snow.data1 <- as.rds.data.frame(HSU.wave1, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")

snow.data2 <- as.rds.data.frame(HSU.wave2, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")

snow.data3 <- as.rds.data.frame(HIV.wave1, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")

snow.data4 <- as.rds.data.frame(HIV.wave2, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")

snow.data5 <- as.rds.data.frame(reRDS, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")

pdf("Resampling_Slovakia.pdf")
reingold.tilford.plot(snow.data1, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia (HSU.Wave1)")
reingold.tilford.plot(snow.data2, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia (HSU.Wave2)")
reingold.tilford.plot(snow.data3, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia (HIV.Wave1)")
reingold.tilford.plot(snow.data4, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia (HIV.Wave2)")

reingold.tilford.plot(snow.data5, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia (RDS)")
reingold.tilford.plot(RDS.SK, vertex.color = "HIV_testing_history", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia")

dev.off()

library(RColorBrewer)
makeboxplot.byway <- function(rep.data, mymain, est, n = 5){
  x0s <- 1:n - 0.4
  x1s <- 1:n + 0.4
  mycol <- brewer.pal(n = n, name = "RdBu")
  
  boxplot(rep.data, col = mycol, xaxt = "n", main = mymain, cex.main = 2, cex.axis = 1.5)
  #mtext("Seeds from", at=0, line = 0.5, side = 1, cex.axis = 1.5)
  #mtext("# of Wave ", at=0, line = 1.5, side = 1, cex.axis = 1.5)
  
  axis(side=1, at=1:4, mgp=c(0,0.5,0), labels = c("HSU", "HSU", "HIVtest", "HIVtest"), cex.axis = 1.5)
  axis(side=1, at=1:4, mgp=c(0,1.5,0), labels = c("Wave1", "Wave2", "Wave1", "Wave2"), cex.axis = 1.5)
  axis(side=1, at=5, mgp=c(0,1,0), labels = c("RDS"), cex.axis = 1.5)
  
  #axis(1, at=1:5,labels=c("HSU.wave1","HSU.wave2", "HIV.wave1", "HIV.wave2", "RDS"), cex.axis = 1)
  segments(x0 = x0s, x1 = x1s, y0 = est, col = "red", lwd = 2)
  legend("topleft", lwd = 2, lty = 1, col = "red", legend = "original RDS data", cex = 1.5)
}


makeboxplot.pdf <- function(n.var, var.name, n.est, pdf.name){
  IT.var <- data.frame(IT.HSU.wave1.rep[,n.var],IT.HSU.wave2.rep[,n.var],IT.HIV.wave1.rep[,n.var],IT.HIV.wave2.rep[,n.var], IT.RDS.rep[,n.var])
  LT.var <- data.frame(LT.HSU.wave1.rep[,n.var],LT.HSU.wave2.rep[,n.var],LT.HIV.wave1.rep[,n.var],LT.HIV.wave2.rep[,n.var], LT.RDS.rep[,n.var])
  RO.var <- data.frame(RO.HSU.wave1.rep[,n.var],RO.HSU.wave2.rep[,n.var],RO.HIV.wave1.rep[,n.var],RO.HIV.wave2.rep[,n.var], RO.RDS.rep[,n.var])
  SK.var <- data.frame(SK.HSU.wave1.rep[,n.var],SK.HSU.wave2.rep[,n.var],SK.HIV.wave1.rep[,n.var],SK.HIV.wave2.rep[,n.var], SK.RDS.rep[,n.var])
  
  pdf(pdf.name)
  makeboxplot.byway(IT.var, mymain = paste("Italy", " (", var.name, ")", sep = ""), est = est.res[1,n.est])
  makeboxplot.byway(LT.var, mymain = paste("Lithuania", " (", var.name, ")", sep = ""), est = est.res[2,n.est])
  makeboxplot.byway(RO.var, mymain = paste("Romania", " (", var.name, ")", sep = ""), est = est.res[3,n.est])
  makeboxplot.byway(SK.var, mymain = paste("Slovakia", " (", var.name, ")", sep = ""), est = est.res[4,n.est])
  dev.off()
}

head(IT.HSU.wave1.rep)
head(IT.RDS.rep1)
makeboxplot.pdf(n.var = 2, var.name = "Age", n.est = 3, pdf.name = "snowball.rep.age.pdf")
makeboxplot.pdf(n.var = 3, var.name = "HSU", n.est = 4, pdf.name = "snowball.rep.HSU.pdf")
makeboxplot.pdf(n.var = 4, var.name = "HIV", n.est = 5, pdf.name = "snowball.rep.HIV.pdf")
makeboxplot.pdf(n.var = 6, var.name = "HIV test history", n.est = 6, pdf.name = "snowball.rep.HIVtest.pdf")
makeboxplot.pdf(n.var = 7, var.name = "ART", n.est = 7, pdf.name = "snowball.rep.ART.pdf")
makeboxplot.pdf(n.var = 8, var.name = "Homosexual", n.est = 8, pdf.name = "snowball.rep.homo.pdf")

###RDSII
makeboxplot.pdf(n.var = 31, var.name = "Age", n.est = 31, pdf.name = "snowball.rep.age.RDSII.pdf")
makeboxplot.pdf(n.var = 9, var.name = "HSU", n.est = 9, pdf.name = "snowball.rep.HSU.RDSII.pdf")
makeboxplot.pdf(n.var = 12, var.name = "HIV", n.est = 12, pdf.name = "snowball.rep.HIV.RDSII.pdf")
makeboxplot.pdf(n.var = 15, var.name = "HIV test history", n.est = 15, pdf.name = "snowball.rep.HIVtest.RDSII.pdf")
makeboxplot.pdf(n.var = 18, var.name = "ART", n.est = 18, pdf.name = "snowball.rep.ART.RDSII.pdf")
makeboxplot.pdf(n.var = 28, var.name = "Homosexual", n.est = 28, pdf.name = "snowball.rep.homo.RDSII.pdf")










makeboxplot.byway(IT.var, mymain = "Italy (Age)", est = est.res[1,3])
makeboxplot.byway(LT.HSU, mymain = "Lithuania (Age)", est = est.res[2,3])
makeboxplot.byway(RO.HSU, mymain = "Romania (Age)", est = est.res[3,3])
makeboxplot.byway(SK.HSU, mymain = "Romania (Age)", est = est.res[4,3])






wave1.IT <- snowball.wave1(RDS.data = RDS.IT, sample.size = dim(RDS.data)[1]
                         , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HSU")
wave1.IT1 <- as.rds.data.frame(wave1.IT, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")
wave1.LT <- snowball.wave1(RDS.data = RDS.LT, sample.size = dim(RDS.data)[1]
                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HSU")
wave1.LT1 <- as.rds.data.frame(wave1.LT, id = "RDS.new.ID"
                               , recruiter.id = "RDS.new.recid"
                               , network.size = "degree")
wave1.RO <- snowball.wave1(RDS.data = RDS.RO, sample.size = dim(RDS.data)[1]
                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HSU")
wave1.RO1 <- as.rds.data.frame(wave1.RO, id = "RDS.new.ID"
                               , recruiter.id = "RDS.new.recid"
                               , network.size = "degree")
wave1.SK <- snowball.wave1(RDS.data = RDS.SK, sample.size = dim(RDS.data)[1]
                           , number.of.coupons = 3, sample.with.replacement = TRUE, seed.var = "HSU")
wave1.SK1 <- as.rds.data.frame(wave1.SK, id = "RDS.new.ID"
                               , recruiter.id = "RDS.new.recid"
                               , network.size = "degree")

pdf("wave1.hsu.pdf")
reingold.tilford.plot(wave1.IT1, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Italy")
reingold.tilford.plot(wave1.LT1, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Lithuania")
reingold.tilford.plot(wave1.RO1, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Romania")
reingold.tilford.plot(wave1.SK1, vertex.color = "HSU", vertex.size = "HIV"
                      , vertex.size.range = c(2,5), main = "Slovakia")
dev.off()

}

