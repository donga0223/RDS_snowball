library(RColorBrewer)

change.colnames <- function(rep.data){
  colnames(rep.data)[27:50] <- paste("RDS", colnames(rep.data)[27:50], sep=".")
  return(rep.data)
}
IT.HIV.wave1.rep <- change.colnames(IT.HIV.wave1.rep)
LT.HIV.wave1.rep <- change.colnames(LT.HIV.wave1.rep)
RO.HIV.wave1.rep <- change.colnames(RO.HIV.wave1.rep)
SK.HIV.wave1.rep <- change.colnames(SK.HIV.wave1.rep)
IT.HIV.wave2.rep <- change.colnames(IT.HIV.wave2.rep)
LT.HIV.wave2.rep <- change.colnames(LT.HIV.wave2.rep)
RO.HIV.wave2.rep <- change.colnames(RO.HIV.wave2.rep)
SK.HIV.wave2.rep <- change.colnames(SK.HIV.wave2.rep)

IT.HSU.wave1.rep <- change.colnames(IT.HSU.wave1.rep)
LT.HSU.wave1.rep <- change.colnames(LT.HSU.wave1.rep)
RO.HSU.wave1.rep <- change.colnames(RO.HSU.wave1.rep)
SK.HSU.wave1.rep <- change.colnames(SK.HSU.wave1.rep)
IT.HSU.wave2.rep <- change.colnames(IT.HSU.wave2.rep)
LT.HSU.wave2.rep <- change.colnames(LT.HSU.wave2.rep)
RO.HSU.wave2.rep <- change.colnames(RO.HSU.wave2.rep)
SK.HSU.wave2.rep <- change.colnames(SK.HSU.wave2.rep)

IT.RDS.rep <- change.colnames(IT.RDS.rep)
LT.RDS.rep <- change.colnames(LT.RDS.rep)
RO.RDS.rep <- change.colnames(RO.RDS.rep)
SK.RDS.rep <- change.colnames(SK.RDS.rep)


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


makeboxplot.pdf <- function(var.name, main.name, pdf.name){
  IT.var <- data.frame(IT.HSU.wave1.rep[[var.name]],IT.HSU.wave2.rep[[var.name]],IT.HIV.wave1.rep[[var.name]],IT.HIV.wave2.rep[[var.name]], IT.RDS.rep[[var.name]])
  LT.var <- data.frame(LT.HSU.wave1.rep[[var.name]],LT.HSU.wave2.rep[[var.name]],LT.HIV.wave1.rep[[var.name]],LT.HIV.wave2.rep[[var.name]], LT.RDS.rep[[var.name]])
  RO.var <- data.frame(RO.HSU.wave1.rep[[var.name]],RO.HSU.wave2.rep[[var.name]],RO.HIV.wave1.rep[[var.name]],RO.HIV.wave2.rep[[var.name]], RO.RDS.rep[[var.name]])
  SK.var <- data.frame(SK.HSU.wave1.rep[[var.name]],SK.HSU.wave2.rep[[var.name]],SK.HIV.wave1.rep[[var.name]],SK.HIV.wave2.rep[[var.name]], SK.RDS.rep[[var.name]])
  
  pdf(pdf.name)
  
  makeboxplot.byway(IT.var, mymain = paste("Italy", " (", main.name, ")", sep = ""), est = est.IT$res[[var.name]])
  makeboxplot.byway(LT.var, mymain = paste("Lithuania", " (", main.name, ")", sep = ""), est = est.LT$res[[var.name]])
  makeboxplot.byway(RO.var, mymain = paste("Romania", " (", main.name, ")", sep = ""), est = est.RO$res[[var.name]])
  makeboxplot.byway(SK.var, mymain = paste("Slovakia", " (", main.name, ")", sep = ""), est = est.SK$res[[var.name]])
  
  if(var.name == ("HIV testing history")){
    var.name <-  "HIV_testing_history"
  }else if(var.name == ("othertesting history")){
    var.name <-  "othertesting_history"
  }
  
  RDS.name <- paste("RDS.", var.name, sep="")
  
  
  IT.RDS.est <- data.frame(IT.HSU.wave1.rep[[RDS.name]],IT.HSU.wave2.rep[[RDS.name]],IT.HIV.wave1.rep[[RDS.name]],IT.HIV.wave2.rep[[RDS.name]], IT.RDS.rep[[RDS.name]])
  LT.RDS.est <- data.frame(LT.HSU.wave1.rep[[RDS.name]],LT.HSU.wave2.rep[[RDS.name]],LT.HIV.wave1.rep[[RDS.name]],LT.HIV.wave2.rep[[RDS.name]], LT.RDS.rep[[RDS.name]])
  RO.RDS.est <- data.frame(RO.HSU.wave1.rep[[RDS.name]],RO.HSU.wave2.rep[[RDS.name]],RO.HIV.wave1.rep[[RDS.name]],RO.HIV.wave2.rep[[RDS.name]], RO.RDS.rep[[RDS.name]])
  SK.RDS.est <- data.frame(SK.HSU.wave1.rep[[RDS.name]],SK.HSU.wave2.rep[[RDS.name]],SK.HIV.wave1.rep[[RDS.name]],SK.HIV.wave2.rep[[RDS.name]], SK.RDS.rep[[RDS.name]])
  
  makeboxplot.byway(IT.RDS.est, mymain = paste("Italy", " (", main.name, ")", sep = ""), est = est.IT$res.RDS[[var.name]][1])
  makeboxplot.byway(LT.RDS.est, mymain = paste("Lithuania", " (", main.name, ")", sep = ""), est = est.LT$res.RDS[[var.name]][1])
  makeboxplot.byway(RO.RDS.est, mymain = paste("Romania", " (", main.name, ")", sep = ""), est = est.RO$res.RDS[[var.name]][1])
  makeboxplot.byway(SK.RDS.est, mymain = paste("Slovakia", " (", main.name, ")", sep = ""), est = est.SK$res.RDS[[var.name]][1])
  dev.off()
}


makeboxplot.pdf(var.name = "age", pdf.name = "snowball.rep.age.pdf", main.name = "Age")
makeboxplot.pdf(var.name = "HIV", pdf.name = "snowball.rep.HIV.pdf", main.name = "HIV")
makeboxplot.pdf(var.name = "HSU", pdf.name = "snowball.rep.HSU.pdf", main.name = "HSU")
makeboxplot.pdf(var.name = "HIV testing history", pdf.name = "snowball.rep.HIVtest.pdf", main.name = "HIV testing history")
makeboxplot.pdf(var.name = "ART", pdf.name = "snowball.rep.ART.pdf", main.name = "ART coverage")
makeboxplot.pdf(var.name = "homosexual", pdf.name = "snowball.rep.homosexual.pdf", main.name = "Homosexual")
makeboxplot.pdf(var.name = "othertesting history", pdf.name = "snowball.rep.othertest.pdf", main.name = "Other testing history")
makeboxplot.pdf(var.name = "howmany_nonsteady_malepartners", pdf.name = "snowball.rep.malepartners.pdf", main.name = "# of Nonsteady Male Partners")
makeboxplot.pdf(var.name = "howmany_nonsteady_malepartners_unprotected", pdf.name = "snowball.rep.malepartners_unprotected.pdf", main.name = "# of Unprotected NSMP")
makeboxplot.pdf(var.name = "injected_drug", pdf.name = "snowball.rep.injected_drug.pdf", main.name = "Injected Drug")
makeboxplot.pdf(var.name = "o1.employed", pdf.name = "snowball.rep.employed.pdf", main.name = "Employed")

###RDSII
makeboxplot.pdf(var.name = "RDS.age", pdf.name = "snowball.rep.RDS.age.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.HIV", pdf.name = "snowball.rep.RDS.HIV.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.HSU", pdf.name = "snowball.rep.RDS.HSU.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.HIV_testing_history", pdf.name = "snowball.rep.RDS.HIVtest.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.ART_coverage", pdf.name = "snowball.rep.RDS.ART.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.homosexual", pdf.name = "snowball.rep.RDS.homosexual.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.other_testing_history", pdf.name = "snowball.rep.RDS.othertest.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.howmany_nonsteady_malepartners", pdf.name = "snowball.rep.RDS.malepartners.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.howmany_nonsteady_malepartners_unprotected", pdf.name = "snowball.rep.RDS.malepartners_unprotected.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.injected_drug", pdf.name = "snowball.rep.RDS.injected_drug.pdf", est.RDS = TRUE)
makeboxplot.pdf(var.name = "RDS.o1.employed", pdf.name = "snowball.rep.RDS.employed.pdf", est.RDS = TRUE)






#############################################################
### Reachable / Not reachable
#############################################################
### reachable sample wave 1

est.samplemean <- function(data){
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
  
  return(res)
}


wave1.reachable <- function(RDS.data, seed.var, t.test.var){
  possible <- union(union(union(union(RDS.data[which(RDS.data[[seed.var]] == 1),"RDS_CODE"],RDS.data[which(RDS.data[[seed.var]]==1),"C1"]),RDS.data[which(RDS.data[[seed.var]]==1),"C2"]),RDS.data[which(RDS.data[[seed.var]]==1),"C3"]), RDS.data[which(RDS.data[[seed.var]]==1),"rec.id"])
  notpossible <- setdiff(RDS.data$RDS_CODE, possible)
  
  possible.RDS <- RDS.data[na.omit(match(possible, RDS.data$RDS_CODE)),]
  notpossible.RDS <- RDS.data[na.omit(match(notpossible, RDS.data$RDS_CODE)),]
  
  possible.RDS$employed <- ifelse(possible.RDS$current_occupation == 1, 1, 0)
  notpossible.RDS$employed <- ifelse(notpossible.RDS$current_occupation == 1, 1, 0)
  
  possiblemean <- est.samplemean(possible.RDS)
  notpossiblemean <- est.samplemean(notpossible.RDS)
  
  test.res <- c()
  for(i in 1:length(t.test.var)){
    test.res[i] <- t.test(possible.RDS[[t.test.var[i]]], notpossible.RDS[[t.test.var[i]]])$p.value
  }
  
  return(list(possiblemean = possiblemean, notpossiblemean = notpossiblemean, test.res = test.res))
}

t.test.var = c("other_testing_history", "howmany_nonsteady_malepartners", "howmany_nonsteady_malepartners_unprotected", "injected_drug", "employed")
reach.IT.HSU.wave1 <- wave1.reachable(RDS.IT, seed.var = "HSU", t.test.var = t.test.var)
reach.LT.HSU.wave1 <- wave1.reachable(RDS.LT, seed.var = "HSU", t.test.var = t.test.var)
reach.RO.HSU.wave1 <- wave1.reachable(RDS.RO, seed.var = "HSU", t.test.var = t.test.var)
reach.SK.HSU.wave1 <- wave1.reachable(RDS.SK, seed.var = "HSU", t.test.var = t.test.var)

reach.IT.HIV.wave1 <- wave1.reachable(RDS.IT, seed.var = "HIV_testing_history", t.test.var = t.test.var)
reach.LT.HIV.wave1 <- wave1.reachable(RDS.LT, seed.var = "HIV_testing_history", t.test.var = t.test.var)
reach.RO.HIV.wave1 <- wave1.reachable(RDS.RO, seed.var = "HIV_testing_history", t.test.var = t.test.var)
reach.SK.HIV.wave1 <- wave1.reachable(RDS.SK, seed.var = "HIV_testing_history", t.test.var = t.test.var)



### reachable sample wave 2

wave2.reachable <- function(RDS.data, seed.var, t.test.var){
  possible1 <- union(union(union(union(RDS.data[which(RDS.data[[seed.var]] == 1),"RDS_CODE"],RDS.data[which(RDS.data[[seed.var]]==1),"C1"]),RDS.data[which(RDS.data[[seed.var]]==1),"C2"]),RDS.data[which(RDS.data[[seed.var]]==1),"C3"]), RDS.data[which(RDS.data[[seed.var]]==1),"rec.id"])
  aa <- union(possible1, na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "C1"]))
  bb <- union(aa,na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "C2"]))
  cc <- union(bb,na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "C3"]))
  possible <- union(cc,na.omit(RDS.data[na.omit(match(RDS.data$RDS_CODE, possible1)), "rec.id"]))
  notpossible <- setdiff(RDS.data$RDS_CODE, possible)
  
  possible.RDS <- RDS.data[na.omit(match(possible, RDS.data$RDS_CODE)),]
  notpossible.RDS <- RDS.data[na.omit(match(notpossible, RDS.data$RDS_CODE)),]
  
  possible.RDS$employed <- ifelse(possible.RDS$current_occupation == 1, 1, 0)
  notpossible.RDS$employed <- ifelse(notpossible.RDS$current_occupation == 1, 1, 0)
  
  possiblemean <- est.samplemean(possible.RDS)
  notpossiblemean <- est.samplemean(notpossible.RDS)
  
  test.res <- c()
  for(i in 1:length(t.test.var)){
    test.res[i] <- t.test(possible.RDS[[t.test.var[i]]], notpossible.RDS[[t.test.var[i]]])$p.value
  }
  
  return(list(possiblemean = possiblemean, notpossiblemean = notpossiblemean, test.res = test.res))
}

t.test.var = c("other_testing_history", "howmany_nonsteady_malepartners", "howmany_nonsteady_malepartners_unprotected", "injected_drug", "employed")
reach.IT.HSU.wave2 <- wave2.reachable(RDS.IT, seed.var = "HSU", t.test.var = t.test.var)
reach.LT.HSU.wave2 <- wave2.reachable(RDS.LT, seed.var = "HSU", t.test.var = t.test.var)
reach.RO.HSU.wave2 <- wave2.reachable(RDS.RO, seed.var = "HSU", t.test.var = t.test.var)
reach.SK.HSU.wave2 <- wave2.reachable(RDS.SK, seed.var = "HSU", t.test.var = t.test.var)

reach.IT.HIV.wave2 <- wave2.reachable(RDS.IT, seed.var = "HIV_testing_history", t.test.var = t.test.var)
reach.LT.HIV.wave2 <- wave2.reachable(RDS.LT, seed.var = "HIV_testing_history", t.test.var = t.test.var)
reach.RO.HIV.wave2 <- wave2.reachable(RDS.RO, seed.var = "HIV_testing_history", t.test.var = t.test.var)
reach.SK.HIV.wave2 <- wave2.reachable(RDS.SK, seed.var = "HIV_testing_history", t.test.var = t.test.var)


#########################################################
#### reachable by variable

reach_byvar <- function(varname, nn, countryname){
  HSU.wave1 <- get(paste("reach.", countryname, ".HSU.wave1", sep = ""))
  HSU.wave2 <- get(paste("reach.", countryname, ".HSU.wave2", sep = ""))
  HIV.wave1 <- get(paste("reach.", countryname, ".HIV.wave1", sep = ""))
  HIV.wave2 <- get(paste("reach.", countryname, ".HIV.wave2", sep = ""))
  print(round(c(HSU.wave1$possiblemean[[varname]], HSU.wave2$possiblemean[[varname]]
          , HIV.wave1$possiblemean[[varname]],  HIV.wave2$possiblemean[[varname]]),4))
  print(round(c( HSU.wave1$notpossiblemean[[varname]],  HSU.wave2$notpossiblemean[[varname]]
                , HIV.wave1$notpossiblemean[[varname]],  HIV.wave2$notpossiblemean[[varname]]),4))
  print(c( HSU.wave1$test.res[[nn]]<0.05,  HSU.wave2$test.res[[nn]]<0.05
                , HIV.wave1$test.res[[nn]]<0.05,  HIV.wave2$test.res[[nn]]<0.05))
  print(c( HSU.wave1$test.res[[nn]]<0.01,  HSU.wave2$test.res[[nn]]<0.01
          , HIV.wave1$test.res[[nn]]<0.01,  HIV.wave2$test.res[[nn]]<0.01))
  print(c( HSU.wave1$test.res[[nn]]<0.001,  HSU.wave2$test.res[[nn]]<0.001
          , HIV.wave1$test.res[[nn]]<0.001,  HIV.wave2$test.res[[nn]]<0.001))
}

reach_byvar("othertesting history", 1, "SK")
reach_byvar("howmany_nonsteady_malepartners", 2, "SK")
reach_byvar("howmany_nonsteady_malepartners_unprotected", 3, "SK")
reach_byvar("injected_drug", 4, "SK")
reach_byvar("o1.employed", 5, "IT")


