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


RDS.rep <- function(RDS.data, setseed = 1001, iter = 100, snowball.sampleway){
  set.seed(setseed)
  snow.data <- snowball.sampleway(RDS.data)
  snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                , recruiter.id = "RDS.new.recid"
                                , network.size = "degree")
  
  HSU <- data.frame()
  aa <- sialon.est(snow.data, snow.RDS)
  HSU <- as.data.frame(cbind(aa$res, aa$res.RDS[1,]))
  
  
  for(p in 2:iter){
    set.seed(setseed+p)
    snow.data <- snowball.sampleway(RDS.data)
    snow.RDS <- as.rds.data.frame(snow.data, id = "RDS.new.ID"
                                  , recruiter.id = "RDS.new.recid"
                                  , network.size = "degree")
    aa <- sialon.est(snow.data, snow.RDS)
    HSU <- rbind(HSU, as.data.frame(cbind(aa$res, aa$res.RDS[1,])))
    print(p)
  }
  return(HSU) 
}


