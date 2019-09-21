#########################################################################################
### Code for Building Logistic regression                                             ###
### Using R function glm()                                                            ###
#########################################################################################

#########################################################################################
### Functions                                                                         ###
###                                                                                   ###
### Dummy: generate dummy variable; first level is the reference level                ###
### Impute: impute missing value                                                      ###
###         Numeric Variable:                                                         ###
###          <= 2 unique values, assign 9999                                          ###
###          if integer, assign mode; other wise assign median                        ###           
###         Categorical Variable:                                                     ###
###         <= 2 unique values, assign "9999"; otherwise assign mode                  ###
### CollapseFactor: group levels of categorical variable; the grouped level is called ###
###                 Other                                                             ###
### getmode:                                                                          ###
### DecileRpt: Generate report                                                        ###
### OverSample: Perform oversample                                                    ###
### BVarSelect: Bootstrapping variable selection                                      ###  
### CompleteSeparate: remove quasi-complete separation variable                       ###
### BuildGlm: Main function to build model                                            ###
######################################################################################### 

#########################################################################################
Dummy <- function(df) {
  for (i in 1:ncol(df)) {   
    x <- df[,i] 
    varname=names(df)[i]    
    x.dummy <- as.data.frame(model.matrix( ~ x -1))
    x.dummy <- x.dummy[-1]    ###### first level is the reference level 
    for (j in 1:length(names(x.dummy))) {
       names(x.dummy)[j] <- gsub("[x]",varname,names(x.dummy)[j])
    }
      
    if (i==1) {
       dummyC <- x.dummy
    } else {
       dummyC <- cbind(dummyC,x.dummy)       
    }
  }  
  return(dummyC)
}
###########################################################################################
Impute <- function(dsn) {
  Names <- names(dsn)
  imputed_value <- data.frame("var"=Names, "Imputed_Value"= rep(NA,length(Names)), "Var_Type"=rep("NA",length(Names)), stringsAsFactors=FALSE)

  for (var in 1:ncol(dsn)) {
   if(length(which(is.na(dsn[,var])))== 0) {
    imputed_value[var,3] <- "NoMissing"
   } else {
     if (is.numeric(dsn[, var])) {
       if (length(unique(dsn[,var])) <= 2){
         imputed <- 9999
         dsn[is.na(dsn[,var]),var] <- imputed    
       } else {
         if (class(dsn[,var])=="integer") {
           imputed <- getmode(dsn[,var])
           dsn[is.na(dsn[,var]),var] <- imputed
         } else {
           imputed <- median(dsn[,var], na.rm = TRUE)
           dsn[is.na(dsn[,var]),var] <- imputed
         } 
       }
       imputed_value[var,2] <- imputed
       imputed_value[var,3] <- "numeric type"
     }
     if (is.character(dsn[, var])) {
       if (length(unique(dsn[,var])) <=2){
         imputed <- "9999"
         dsn[is.na(dsn[,var]),var] <- imputed
       } else {
         imputed <- getmode(dsn[,var])
         dsn[is.na(dsn[,var]),var] <- imputed
       }
       imputed_value[var,2] <- imputed 
       imputed_value[var,3] <- "character_type"
     }
    }
  }
  return(list(dsn, imputed_value))
}
######################################################################################
getmode <- function(v) {
  v <- v[!is.na(v)]
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#######################################################################################
DecileRpt <- function(x,group) {  
  score <- x[,"score"]
  Y <- x[,2]  
  slice <- cut2(score,g=group,digits=9)
  min_prob <- tapply(score,slice,min)
  pred_prob <- tapply(score,slice, mean)
  circ <- table(slice)
  est_target <- as.integer(pred_prob * circ)
  actual_target <- tapply(Y,slice, sum)
  actual_prob <- (actual_target/circ)
  slice_rpt <- cbind(min_prob, circ,actual_target,actual_prob,pred_prob,est_target)
  slice_rpt <- slice_rpt[order(-slice_rpt[,1]),]
  decile <- seq(1,dim(slice_rpt)[1])
  slice_rpt <- cbind(decile,slice_rpt) 
  rownames(slice_rpt) <- decile
  return(slice_rpt)
}
#########################################################################################
BVarSelect <- function(dsn, boot.num=5,keep.num=2, sig_level=0.1) {
  for (run in 1:boot.num) {
    size <- 0.5 + runif(1)/2
    smpl.size <- ceiling(nrow(dsn) * size)
    smpl.index <- sample(1:nrow(dsn), smpl.size, replace=FALSE)
    dsn.sample <- dsn[smpl.index,]
    fit <- glm(dsn.sample[,1] ~ ., family=binomial(link="logit"),data=dsn.sample[,-1])
    fit.m <- summary(fit)$coeff
    if (run==1) {
      fit.3 <- fit.m  
    } else {
      fit.3 <- rbind(fit.3, fit.m)
    }
  } 
  ind <- which(rownames(fit.3)=="(Intercept)")
  fit.3 <- fit.3[-ind,]
  rm(ind)
  keep.var <- rownames(fit.3[fit.3[,4] <= sig_level,])
  keep.final <- names(table(keep.var)[table(keep.var) >= keep.num])
  keep_final <- unique(keep.final)
  return(keep_final)
}
##############################Over Sample########################################################
OverSample <- function(df, resp_rate=0.2, dep_var=dep.var) {
  p <- mean(df[,dep_var],na.rm=TRUE)
  old.odds <- log(p/(1-p))
  good.index <- which(df[dep_var]==1)
  good <- df[good.index,]
  bad <- df[-good.index,]
  good.obs=nrow(good)
  bad.obs <- nrow(bad)
  big_good <- good
  if (p <= resp_rate) { 
    required_good_size <- ceiling(bad.obs * (resp_rate/(1-resp_rate)))
    replica <- ceiling(required_good_size/good.obs)
    if (replica == 1) {
      big_good <- good
    } else {
      for (i in 1:(replica - 1)) {
        big_good <- rbind(big_good,good)
      }
    }  
    smpl.size <- required_good_size
  } else {
    smpl.size <- good.obs
  }
  big_good.obs <- nrow(big_good)
  smpl.index <- sample(1:big_good.obs, min(smpl.size,big_good.obs), replace=FALSE)
  good <- big_good[smpl.index,]
  df <- rbind(good,bad)
  new_p <- mean(df[,dep_var], na.rm=TRUE)
  new.odds <- log(new_p/(1 - new_p))
  adjustment.value <- old.odds - new.odds
  rm(p, good.index, good, bad, smpl.size, smpl.index, bad.obs,good.obs,big_good.obs,
     required_good_size,replica,i,new_p)
  return(list(df,adjustment.value,old.odds,new.odds))
}
##################################################################################################
CollapseFactor <- function(x,namex,p) {
     t <- table(x)
     prop <- t/sum(t)
     levt <- cbind(names(t), names(t)) 
     ind1<- which(prop < p)
     
     if (length(ind1) > 0) {
       prop2 <- prop[-ind1]
       ind1 <- names(prop)[ind1]
       ind2 <- which.min(prop2)
       ind2 <- names(prop2)[ind2]
       ind <- c(ind2,ind1)
       ind <- which(levt[,1] %in% ind)
       rm(ind1,ind2)

       levt[ind, 2] <- "Other"
       collapseLevel <- paste(names(t)[ind],collapse=",")
       collapseVar <- namex
       lv <- levels(as.factor(x))
       if(sum(sort(lv) != sort(levt[, 1]))>0) {
        stop ("The names from substitution table does not match given level names")
       }
       res <- rep(NA, length(x))
       for(j in lv) {
         res[x==j] <- as.character(levt[levt[, 1]==j, 2])
       }
     } else {
       collapseLevel <- "NO_COLLAPSE"
       collapseVar <- namex
       res <- x 
     }
     out_list <- list(factor(res), collapseVar, collapseLevel)
     return(out_list)
}
############################################################################################
CompleteSeparate <- function(df, catvar, targetvar, bound=0.9) {
  drop <- NULL
  df <- df[,catvar,drop=FALSE]
  for (i in 1:length(catvar)) {
     t <- as.data.frame(as.table(tapply(target,df[,i],sum,na.rm=TRUE)))
     t$pct <- t[,2]/sum(t[,2])
     if (any(t$pct >= bound)) {
        droptemp <- catvar[i]
     } else {
        droptemp <- NULL
     }
     drop <- c(drop,droptemp)
     rm(droptemp)
  }
  return(drop)
}  
############################################################################################
############################################################################################
BuildGlm <- function(ds, dep.var, num_treat_as_cat="Y", cutoff_unique_value=5, oversample="N",prop_cutoff=0.02,
                     desired_resp_rate=0.2, split_size=0.67,p_value=0.05,vif_cutoff=4, bootNum=5,keepNum=2,
                     lower_bound=0.05,upper_bound=0.95,SeparateBound=0.9) {

  set.seed(as.integer(format(Sys.Date(),"%Y%m%d"))) 
  Y  <- ds[dep.var]
  ds <- ds[,-1]  ##### exclude dep.var from dsn #####
############################################################################################
  CountNA_G<-sapply(ds,function(y)length(which(is.na(y)==T)))
  Count_Miss<- data.frame(model_var=colnames(ds), Count=CountNA_G, MissPCT = CountNA_G/nrow(ds)) 
#############################################################################################
### Missing: high (>=0.75); midium (0.25 t0 0.7499);low (< 0.25) ############################
  MissAll <- names(ds[CountNA_G/nrow(ds) == 1.0 ] ) 
  MissHigh <- names(ds[(CountNA_G/nrow(ds) >= .75 ) & (CountNA_G/nrow(ds) < 1 )] ) 
  MissMedium <- names(ds[(CountNA_G/nrow(ds) >= .25 )&  (CountNA_G/nrow(ds) < .75 )] ) 
  MissLow <- names(ds[(CountNA_G/nrow(ds)< .25 )&  (CountNA_G/nrow(ds) >0 )] ) 

  NoMiss <-  names(ds[CountNA_G/nrow(ds) == 0.0] ) 

  Missing <- with (Count_Miss, ifelse(MissPCT ==1, "MissALL", ifelse(MissPCT >= 0.75,"MissHigh", 
                   ifelse(MissPCT >= 0.25, "MissMedium", ifelse(MissPCT > 0, "MissLow", "NoMiss")))))

  Count_Miss$Missing_Level <- Missing
############# Save Missing Value Statistics ####################################################
  write.csv(Count_Miss,file="Missing_Var_Statistics.csv")
################################################################################################
  ds <- ds[,!(names(ds) %in% c(MissAll, MissHigh))]
  rm(CountNA_G, Count_Miss, Missing)
################################################################################################
#################### Numeric Variables                          ################################
  num_list<- sapply(ds, is.numeric)
################################################################################################
#################### categorical Variables                      ################################
  factor_var <- colnames(ds)[!num_list]
  if (length(factor_var) > 0) {
    for (i in factor_var) {
      ds[,i] <- as.character(ds[,i])
    }
  }
###############################################################################################
###################### Imputed  Missing #######################################################
  junk <- Impute(ds)
  ds <- junk[[1]]
  imputed_value <- junk[[2]]
  rm(junk)
###############################################################################################
############### Numeric Variables #############################################################
  num_list<- sapply(ds, is.numeric)
  levels_num <- NULL
  dsNumVar <- NULL
  if (length(num_list[num_list==TRUE]) > 0) {
    num_vars<- colnames(ds)[num_list]
    dsNumVar <- ds[num_vars]
    levels_num<-sapply(ds[num_vars], function(y) length(unique(y)))
######### Single value var ####################################################################
    single_value_var <- names(levels_num[levels_num <= 1])
    levels_num <- levels_num[!(names(levels_num) %in% single_value_var)]
###### Binary value Var and exclude biased category var #######################################
    binary_num <- names(levels_num[levels_num == 2]) 
    bDrop <- NULL
    if (length(binary_num) > 0) {
      dst <- ds[binary_num]
      for (i in 1:ncol(dst)) {
        x <- prop.table(table(dst[,i]))
        if (any(x < prop_cutoff)) {
         bDrop <- c(bDrop,names(dst)[i])
        }
      }
      rm(dst)
    }
    levels_num <- levels_num[!(names(levels_num) %in% bDrop)]
    
    if (length(levels_num) > 0) {
      dsNumVar <- dsNumVar[names(levels_num)]
    } else {
      dsNumVar <- NULL
      numVar <- NULL
    }

    if (num_treat_as_cat=="Y" & length(levels_num) > 0 ) {
      cat_var_from_num <- names(levels_num[levels_num <= cutoff_unique_value])
    } else {
      cat_var_from_num <- NULL
    } 
    
    if (nrow(dsNumVar) > 0) {
      dsNumVar <- dsNumVar[!(names(dsNumVar) %in% cat_var_from_num)]
      numVar <- names(dsNumVar)  
    }
   
    rm(bDrop)
    rm(single_value_var,num_vars, binary_num)
  }
#################################################################################################
#### Bound Numeric Variable #####################################################################
  if (nrow(dsNumVar) > 0) {
    pLower <- apply(dsNumVar,2, function(x) quantile(x,probs=lower_bound,na.rm=TRUE)) 
    pUpper <- apply(dsNumVar,2, function(x) quantile(x,probs=upper_bound,na.rm=TRUE))

    for (i in 1:length(names(dsNumVar))) {
      dsNumVar[,i] <- pmax(pmin(dsNumVar[,i],pUpper[i]),pLower[i])
    }

    pBound <- data.frame(lowerBound=pLower,upperBound=pUpper)
    var <- row.names(pBound)
    pBound <- cbind(var,pBound)
    row.names(pBound) <- NULL
  } 
#################################################################################################
############################## Non Numeric Variables             ################################
  factor_var <- colnames(ds)[!num_list]
  factor_var <- unique(c(factor_var,cat_var_from_num))
#################################################################################################
###### remove quasi separation variable                          ################################ 
  sepvar <- CompleteSeparate(ds,factor_var,Y[,1],bound=SeparateBound)
  factor_var <- factor_var[!(factor_var %in% sepvar)]
  rm(sepvar)

##################################################################################################  
  catVar <- NULL
  dsBCharVar <- NULL
  dsMCharVar <- NULL
  if (length(factor_var) > 0) {
    dsCat <- ds[factor_var]
    for (i in factor_var) {
      dsCat[,i] <- as.character(dsCat[,i])
    }
#################################################################################################
    levels_char <- sapply(dsCat, function(y) length(unique(y)))
############### Single value Char Var  ##########################################################
    single_value_char <- names(levels_char[levels_char <= 1])
############## Multi Value Char Var #############################################################
    multi_value_char  <-  names(levels_char[levels_char > 2])
####### Binary Value var #######################################################################
    binary_value_char <- names(levels_char[levels_char == 2]) 
#################################################################################################
######################## Exclude biased Binary category var #####################################
    bDrop <- NULL

    cat_var_from_binary_char <- NULL

    if (length(binary_value_char) > 0) {
      cat_var_from_binary_char <- factor_var[factor_var %in% binary_value_char]
      dst <- dsCat[binary_value_char]
      for (i in 1:ncol(dst)) {
        x <- prop.table(table(dst[,i]))
        if (any(x < prop_cutoff)) {
          bDrop <- c(bDrop,names(dst)[i])
        }
      } 
      cat_var_from_binary_char <- cat_var_from_binary_char[!(cat_var_from_binary_char %in% bDrop)]
      rm(dst,bDrop,i,binary_value_char,x)
    }
    
    if (length(cat_var_from_binary_char) > 0) {
      dsBCharVar <- dsCat[cat_var_from_binary_char]
    } else {
      dsBCharVar <- NULL
    } 
#################################################################################################
############## Multi Value Category Var #########################################################
    cat_var_from_multi_char <- NULL
 
    if (length(multi_value_char) > 0 ) {
      dsCatM <- dsCat[multi_value_char]
      nameChk <- names(dsCatM)
      cat_var <- character(ncol(dsCatM))
      collapse <- character(ncol(dsCatM))
      for (k in 1:ncol(dsCatM)) {
        w <- CollapseFactor(dsCatM[,k],nameChk[k], prop_cutoff)
        dsCatM[,k] <- w[[1]]
        cat_var[k] <- w[[2]]
        collapse[k] <- w[[3]]
      }
      map_cat <- as.data.frame(cbind(cat_var,collapse))
      names(map_cat) <- c("Var", "Levels_Grouped_Together")
      rm(cat_var,collapse,nameChk,w,k)
##############################################################################################
      levels_MValue_char <- sapply(dsCatM, function(y) length(unique(y)))
######### Single value var ####################################################################
      single_value_MVar <- names(levels_MValue_char[levels_MValue_char <= 1])
      cat_var_from_multi_char <- names(levels_MValue_char)[!(names(levels_MValue_char) %in% single_value_MVar)]
      map_cat <- map_cat[map_cat[,1] %in% cat_var_from_multi_char,]

      if (length(cat_var_from_multi_char) > 0 ) {
        dsMCharVar <- dsCatM[cat_var_from_multi_char] 
      } else {
        dsMCharVar <- NULL
      }

      rm(multi_value_char,dsCatM,single_value_MVar,levels_MValue_char) 
    }  
    catVar <- c(cat_var_from_binary_char, cat_var_from_multi_char)
  }
###########################################################################################
################################################################################################
  if (length(nrow(dsBCharVar)) > 0) {
     if (length(nrow(dsMCharVar)) > 0) {
        junk <- cbind(dsBCharVar,dsMCharVar)
     } else {
        junk <- dsBCharVar
     }
  } else {
    if (length(nrow(dsMCharVar)) > 0) {
      junk <- dsMCharVar
    } else {
      junk <- NULL
    }
  }

   if (length(nrow(junk)) > 0) {
     if (length(nrow(dsNumVar)) > 0) {
        ds2 <- cbind(junk,dsNumVar)
     } else {
        ds2 <- junk 
     }
   } else {
     if (length(nrow(dsNumVar)) > 0) {
       ds2 <- dsNumVar
     } else {
       ds2 <- NULL
     }
   }
  if (is.null(ds2)) stop("No Significant Predictors") 
  rm(dsBCharVar,dsMCharVar,dsNumVar,junk)

### Modeling data set ##########################################################################
  ds2 <- cbind(Y,ds2)
  
####################Over Sample   ############################################################### 
  if (oversample=="Y") {
    inputdata <- ds2 
    junk <- OverSample(ds2,resp_rate= desired_resp_rate,dep_var=dep.var)
    ds2=junk[[1]]
    adjustment_value=junk[[2]]
    rm(junk)
  }
#################################################################################################
##################### Variable Selection   ######################################################
  tgtVar <- ds2[1]
  if (length(numVar) >=1) {
    testN <- cbind(ds2[1],ds2[numVar])
  } else {
    testN <- NULL
  }
  if (length(catVar) >=1) {
    testC <- cbind(ds2[1],ds2[catVar])
  } else {
    testC <- NULL
  }
  rm(ds2)

##################### Numeric variable     ######################################################
  if (length(testN) > 2) {
    test <- apply(testN[,-1],2, function(x,y=as.numeric(testN[,1])) {
         if (length(unique(x)) > 1) {
           fit <- lm(y ~ x)
           jk <- summary(fit)$fstatistic
           output <- return(c(summary(fit)$r.squared,jk[1],pf(jk[1],jk[2],jk[3],lower.tail=F)))       
         } else {
           output <- return(c(0.0, 0.0, 10.0))
         }
      }
    ) 

    test <- as.data.frame(t(test))
    names(test) <- c("r.square","f.stat","p.value")
    var <- row.names(test)
    test <- cbind(var,test)
    row.names(test) <- NULL
 
    keep_var_n <- as.character(test[test$p.value < p_value,][,1])
    if (length(keep_var_n) > 0) {
      testN <- testN[,keep_var_n]
    } else {
      testN <- NULL
    }
    rm(test)
  } else {
    if (length(testN) == 2) {
      keep_var_n <- names(testN)[2]
    } else {
      keep_var_n <- NULL
    }
  }  
################## Category var ###############################################################
  if (length(testC) > 2 ) {
    test <- apply(testC[,-1],2, function(x,y=as.numeric(testC[,1])) {
         x <- as.factor(x)
         if (length(levels(x)) > 1) {
           x.dummy <- as.data.frame(model.matrix(~ x -1))
           x.dummy <- x.dummy[-1] 
           fit=lm(y ~ ., data=x.dummy)
           jk <- summary(fit)$fstatistic
           output <- return(c(summary(fit)$r.squared,jk[1],pf(jk[1],jk[2],jk[3],lower.tail=F)))       
         } else {
           output <- return(c(0.0, 0.0, 10.0))
         }
       }
    ) 

    test <- as.data.frame(t(test))
    names(test) <- c("r.square","f.stat","p.value")
    var <- row.names(test)
    test <- cbind(var,test)
    row.names(test) <- NULL

    keep_var_c <- as.character(test[test$p.value < p_value,][,1])
    if (length(keep_var_c) > 0) {
      testC <- testC[keep_var_c]
      for (i in keep_var_c) {
        testC[,i] <- as.factor(testC[,i])
      }
    } else {
      testC <- NULL
    }  
    rm(test)
    
    if (length(testC) > 0) {
      testC <- as.data.frame(Dummy(testC))
    }
  } else {
    if (length(testC) == 2) {
      keep_var_c <- names(testC)[2]
      testC[,2] <- as.factor(testC[,2])
      testC <- testC[2]
      testC <- as.data.frame(Dummy(testC))     
    } else {
      keep_var_c <- NULL
    }    
  }
##########################################################################################
######################## Write imputed_vale, pBound, map_cat #############################
  imputed_value <- imputed_value[(imputed_value[,1] %in% c(keep_var_n,keep_var_c)),]
  write.csv(imputed_value, file="Missing_Value_Imputed.csv",row.names=FALSE)

  pBound <- pBound[(pBound[,1] %in% keep_var_n),]
  write.csv(pBound, file="Bound_of_Num_Var.csv", row.names=FALSE)  
  write.csv(map_cat, file="Collapse_of_Cat_Var.csv",row.names=FALSE)
#############################################################################################
  if (ncol(testC) > 0) {
    ds2 <- cbind(tgtVar,testC)
  } else {
    ds2 <- tgtVar
  }

  if (ncol(testN) > 0) {
    ds2 <- cbind(ds2,testN)
  } else {
    ds2 <- ds2
  }
############ check Multivollinearity ##########################################################
  xvif <- vif(lm(ds2[,1] ~ . , data=ds2[,-1]))
  ind <- which(xvif <= vif_cutoff)
  keep_var <- as.data.frame(xvif[ind])
  keep_name <- row.names(keep_var)
  for (i in 1:length(keep_name)) {
    keep_name[i] <- gsub('\\`','',keep_name[i])
  }
  keep_name <- c(dep.var,keep_name)
  
  if (length(keep_name) <= 1) stop("No Significant Predictors") 

  ds2 <- ds2[,keep_name]
  rm(keep_name)
##############################################################################################
################# Split ######################################################################
  smpl.size <- ceiling(nrow(ds2) * split_size)
  smpl.index <- sample(1:nrow(ds2), smpl.size, replace=FALSE)
  train <- ds2[smpl.index,]
  validate <- ds2[-smpl.index,]
################################################################################################
  keep_var <- BVarSelect(train, boot.num=bootNum, keep.num=keepNum, sig_level=p_value)
  
  if (length(keep_var) <= 1) stop("No Significant Predictors") 
################################################################################################
  for (i in 1:length(keep_var)) {
    keep_var[i] <- gsub('\\`','',keep_var[i])
  }

  train <- train[names(train) %in% c(dep.var,keep_var)]

###### Keep vatiables in training dataset ######################################################
  trainVar <- names(train)[-1]

################################################################################################
  glm_fit <- glm(train[,1] ~ . , family=binomial(link="logit"),data=train[,-1])
  fit <- summary(glm_fit)$coeff
  rownames(fit)[1] <- "intercept"
  model <- fit
  model <- as.data.frame(model)
  rm(fit,sample)
  model_var <- rownames(model)[-1]  
  estimate=model[,1]
  for (i in 1:length(model_var)) {
    model_var[i] <- gsub('\\`','',model_var[i])
  }
##############################################################################################
  score <- predict(glm_fit,train,type="response")
  score_y <- cbind(score,train[,1])
###### Plot ROC Curve ##########################################################################
  actuals <- train[,1][order(score)]
  sens <- (sum(actuals) - cumsum(actuals))/sum(actuals)
  spec <- cumsum(!actuals)/sum(!actuals)
  plot(1 - spec,sens,type="l",col="red", ylab="Sensitivity", xlab="1 - Specificity", main="ROC Curve")
  abline(c(0,0),c(1,1))
  print("Decile Report and auc : Train Data")
  decile_rpt <- DecileRpt(score_y,10)
  auc <- sum(spec * diff(c(0,1-sens)))
  print(decile_rpt)
  print(auc)
############### Validate ###############################################################
  score <- predict(glm_fit,validate,type="response")
  score_y <- cbind(score,validate[,1])
  print("Decile Report: Validate Data")
  decile_rpt <- DecileRpt(score_y,10)
  print(decile_rpt)
##########################################################################################  
  dsnVar <- train[model_var]
  allone <- rep(1,nrow(train))
  dsnVar <- cbind(allone,dsnVar)
  Standardized_Estimate <- (apply(dsnVar,2, sd, na.rm=TRUE)/(pi/sqrt(3)))
  Standardized_Estimate <- Standardized_Estimate * estimate
  Var_Weight <- abs(Standardized_Estimate)/sum(abs(Standardized_Estimate), na.rm=TRUE)
  Var_Weight[1] <- NA
  model_var <- as.data.frame(c("intercept",model_var))
  model <- cbind(model_var,model,Var_Weight)
  names(model)[1] <- "Variable"
  rownames(model) <- NULL
################## Score Input data #######################################################
  if (oversample=="Y") {
    testN <- NULL
    if (length(keep_var_n) > 0) {
      testN <- inputdata[keep_var_n]
    }
    testC <- NULL
    if (length(keep_var_c) > 0) {
      testC <- inputdata[keep_var_c]
      for (i in keep_var_c) {
        testC[,i] <- as.factor(testC[,i])
      }
      testC <- as.data.frame(Dummy(testC))
    }
    if (length(testC) > 0 & length(testN) > 0) {
       inputdata1 <- cbind(testC,testN)
    }

    if (length(testC) > 0 & length(testN) ==0) {
       inputdata1 <- testC
    }
    if (length(testC) == 0 & length(testN) > 0) {
      inputdata1 <- testN
    }

    for (i in 1:length(names(inputdata1))) {
      names(inputdata1)[i] <- gsub('\\`','',names(inputdata1)[i])
    }

    ds4 <- inputdata1[,as.character(model[,1])[-1]]
    rm(inputdata1) 
    intercept <- rep(1, nrow(ds4))
    ds4 <- cbind(intercept, ds4)
    estimate[1] <- estimate[1] + adjustment_value
    logit <- as.matrix(ds4) %*% estimate
    colnames(logit) <- "logit"
    score <- exp(logit)/(1+exp(logit))
    colnames(score) <- "score"
    score_y <- cbind(score,inputdata[,1])
    decile_rpt <- DecileRpt(score_y,10)
    print("Decile Report: After adjusting Over Sample")
    print(decile_rpt)
    model[1,2] <- model[1,2] + adjustment_value
    rm(testN,testC,ds4,score_y,intercept)
  } 
  print("Fitted Model")
  print(model)

}
###########################################################################################
###########################################################################################






#########################################################################################
### Parameters:                                                                       ###
### ds : input data                                                                   ###
### dep.var : target variable                                                         ###
### oversample: "Y": want to oversample; "N": do not want to oversample               ###
### desired_resp_date: if oversample, the desired response rate in the oversample file###
### prop_cutoff: Threshold for dropping a variable having binary value only;          ###
###              if the percentage of any one value  < prop_cutoff, the variable will ###
###              be dropped                                                           ### 
### num_treat_as_cat: "Y": want to treat the numeric variables having value <=        ###  
###                   cutoff_unique_value as categorical variable                     ###
###                   "N": do not want treat any numeric variable as categorical var  ###
### cutoff_unique_value                                                               ###
### split_size:                                                                       ###
### p_value:                                                                          ###
### vif_cutoff: threshold for measuring multillinearity                               ###
### bootNum: numer of bootstrapping                                                   ###
### keepNum: threshold to keep variables at the end of bootstrapping runs             ###
### lower_bound: lower percentile for treating outlier                                ###
### upper_bound: upper percentile for treating outlier                                ###
### SeparateBound: lower bound for removeing quasi-complete separation variable       ###  
#########################################################################################

########### Read in Data and Build the model ##############################################
############################################################################################
options(warn=-1)
require(Hmisc)
require(DAAG)

#setwd("L:/DS_Modeling_Competition/Comp1")
setwd("c:/Users/ychen/Nov7_2016/Yr2018/SBA_Modeling_Competition/Census_Income_Dataset/Classification_1_Binary")

ds <- read.csv("census_train.csv", na.strings=c("",".",NA,NaN,"NA"," "), stringsAsFactors=FALSE)

names(ds) <- tolower(names(ds)) 

target <- ifelse(ds$target_income=="50000+.",1,0)

ds <- ds[,-c(1,43)]

ds <- cbind(target,ds)
dep.var <- "target"

BuildGlm(ds, dep.var, 
         oversample="N", desired_resp_rate=0.2 ,prop_cutoff=0.01,
         num_treat_as_cat="Y", cutoff_unique_value=2, 
         split_size=0.67, p_value=0.1, vif_cutoff=4,
         bootNum=10, keepNum=5,lower_bound=0.01,upper_bound=0.99, SeparateBound=0.9)





