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

########################################################################################
########### as.RDS.data ################################################################
########################################################################################
########################################################################################

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


  
########################################################################################
########### sample and RDS estimations #################################################
########################################################################################
########################################################################################

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

  
########################################################################################
########### Snowball sampling ##########################################################
########################################################################################
########################################################################################

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


snowball.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway, seed.var = seed.var){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data, seed.var = seed.var)
  snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
  
  HSU <- data.frame()
  aa <- sialon.est(snow.data, snow.RDS)
  HSU <- as.data.frame(cbind(aa$res, aa$res.RDS[1,]))
  
  #if(length(which(snow.data1$HIV==1))==0){
  #  a2 <- a1
  #  a2$estimate <- a2$interval[2] <- a2$interval[3] <- 0
  #}else{
  #  a2 <- RDS.II.estimates(rds.data=snow.data1 ,outcome.variable='ART_coverage', subset = HIV == 1)
  #}
  for(i in 2:iter){
    set.seed(setseed+i)
    snow.data <- snowball.sampleway(RDS.data, seed.var = seed.var)
    snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
    aa <- sialon.est(snow.data, snow.RDS)
    HSU <- rbind(HSU, as.data.frame(cbind(aa$res, aa$res.RDS[1,])))
    print(i)
  }
  return(HSU) 
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

IT.RDS.rep <- snowball.rep(RDS.data = RDS.IT, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)
LT.RDS.rep <- snowball.rep(RDS.data = RDS.LT, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)
RO.RDS.rep <- snowball.rep(RDS.data = RDS.RO, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)
SK.RDS.rep <- snowball.rep(RDS.data = RDS.SK, setseed = 1001, iter = 1000, snowball.sampleway = RDS.sample.dongah)

print("6.done")
save.image("sialon_all_rep.RData")




