# ---
# This script is a companion to the main script titled "SDHMmain"
# This script stores functions called by the main script and should not be altered.
# Written by William Wiskes
# Last update 9/6/2022
# ---

#this function will pull in data from a postgres database
library(RPostgreSQL)
queryPostgres <- function(species) {
    source("env.R") #set all the variables needed from a postgres database in an external env file
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)  
    query <- dbSendQuery(con, paste0("SELECT * FROM ",species))
    occs <- dbFetch(query)
    occs[,1:3] %>% drop_na()
}


library(bigrquery)
queryBiobase <- function(species) {
  # Store the project ID
  projectid = "ut-dnr-biobase-dev"
  
  # Set your query
  sql <- paste0("SELECT * FROM `ut-dnr-biobase-dev.biobase.", species, "`")
  
  # Run the query; this returns a bq_table object that you can query further
  tb <- bq_project_query(projectid, sql)
  
  # Store the data in a tibble
  out <-bq_table_download(tb)
  
  out
}

#This is a function to loop through your raster list and punch through your occurrence data
library(rgdal)
extractStack <- function(pointData, rasterList) {
    i = 1
    d = pointData
    len = length(rasterList)
    while (i < len+1){
        #a <- substr(rasterList[i],5,nchar(rasterList[i]))
        a <- rasterList[i]
        b <- paste0("/vsicurl/https://storage.googleapis.com/predictors_public/",a)
        c <- stack(b)
        d <- cbind(d, raster::extract(c, pointData[,c("x","y")]))
        i = i + 1
    }
d
}

#To generate pseudo absences from the extent of the study area
library(sfheaders)
pseudoFunction <- function(pointSF, blob, temp) {
    #intersect points with study area, this is redundant with state data, but may be important for other workflows
    SF <- pointSF %>% st_intersection(blob)
    #make a dataframe
    aea <- sf_to_df(SF)
    #subset to just the lat/lng
    aea <- aea[tail(names(aea), 2)]
    #add sppres column
    aea$sppres <- 1
    #rename the lat/lng columns
    names(aea)[1:2] <- c("aea_xF","aea_yF")

    #set sampling frame
    blob$in_frame <- 1
    frame <-rasterize(blob,temp, "in_frame")

    #Make fnet
    f1 <- sp::coordinates(temp) # get spatial coords 
    f2 <- cellFromXY(temp, f1) # grab cell number for each [X,Y] coordinate value
    fnetDF <- as.data.frame(cbind(f2, f1)) # build dataframe
    names(fnetDF) <- c("FNETID", "cell_aea_x", "cell_aea_y")

    ##### pseudo abs

    tru.xy <- aea[c("aea_xF","aea_yF")]
    p1 <- coordinates(tru.xy)
    fnetid <- cellFromXY(frame, p1)

    FNET <- cbind(fnetid, aea)
    names(FNET)[1] <- "FNETID"

    #pseudo pres
    fnetid <- raster::extract(frame, fnetDF[c("cell_aea_x","cell_aea_y")])

    #merge
    tfnetid <- cbind(fnetDF, fnetid)
    m1 <- merge(fnetDF, tfnetid, by = c("FNETID", "cell_aea_x", "cell_aea_y"), all.y = T)

    indexfnet <- merge(m1, FNET, by = c("FNETID"), all = T)

    names(indexfnet)[4] <- "in.modFR"

    #draw pseudo absence

    p1<- indexfnet
    p2.spp <- subset(p1$FNETID, p1$sppres == 1)
    p2.modFR <- subset(p1$FNETID, p1$in.modFR == 1)

    p3 <- p1[!p1$FNETID %in% p2.spp, ]
    p4 <- p1[!p1$FNETID %in% p2.spp & p1$FNETID %in% p2.modFR, ]

    set.seed(1234)
    srs <- p4[sample(1:nrow(p4), 2 * table(p1$sppres)[[1]], replace = F), ] #make this an input variabe with a default of two

    srs$sppres <- 0
    srs$in.modFR <- 1

    #merge pseudo abs
    PPSA <- merge(FNET, srs, by = c("FNETID", "sppres", "aea_xF", "aea_yF"), all = T)
    PPSA$in.modFR <- NULL


    PPSA$x <- ifelse(PPSA$sppres == 0, PPSA$cell_aea_x, PPSA$aea_xF)


    PPSA$y <- ifelse(PPSA$sppres == 0, PPSA$cell_aea_y, PPSA$aea_yF)

    PPSA
}

#This function removes columns from the data correlated to eachother according to spearmans 
#additional remove or preserve columns are added/removed based on lists within their respective variable names
cutFunction <- function(data, cut, preserve, remove) {
    #remove non-data columns
    if(nchar(remove[1])>0){
    exclude <- names(data) %in% c('FNETID', 'sppres', 'aea_xF','aea_yF', 'cell_aea_x','cell_aea_y','x','y', remove)  
        } else {
    exclude <- names(data) %in% c('FNETID', 'sppres', 'aea_xF','aea_yF', 'cell_aea_x','cell_aea_y','x','y')
        }
    #test correlation 
    df1 <- cor(data[!exclude], use = "pairwise.complete.obs", method = "spearman") 
    #remove na rows
    df2 <- as.data.frame(df1) %>% filter(if_any(everything(), ~!is.na(.)))
    #remove na columns
    all_na <- function(x) any(!is.na(x))
    df3 <- df2 %>% select_if(all_na)
    noNA <- colnames(df3)
    #return a dataframe which all that have a correlation above the abs of cut as TRUE
    df4 <-as.data.frame(subset(df3 > cut | df3 < -cut))
    BadCol <- c()
    #append column names over cut threshold to a list except the first variable correlated to them.
    #only check subsiquent columns if not already in list.
    
    # for (x in colnames(df4)) {
    #   var <- df4 %>%
    #     dplyr::select(x) %>%
    #     filter(get(x) != FALSE)
    #   if (row.names(var) != colnames(var)) {
    #     #print(row.names(var)[2:length(row.names(var))])
    #     if (!(colnames(var) %in% BadCol)){
    #     BadCol <- append(BadCol,row.names(var)[2:length(row.names(var))])
    #     }
    #   }
    #   
    #   #print(get(paste0('df4$', x)))
    # }
    
    for (x in colnames(df4)) {
      var <- df4 %>%
        dplyr::select(x) %>%
        filter(get(x) != FALSE)
      var = var[row.names(var) != colnames(var), , drop = FALSE]
      #print(var)
      #if (row.names(var) != colnames(var)) {
      for (x in row.names(var)) {
        #print(row.names(var)[2:length(row.names(var))])
        if (!(colnames(var) %in% BadCol)){
          #BadCol <- append(BadCol,row.names(var)[2:length(row.names(var))])
          BadCol <- append(BadCol,x)
        }
      }
      
      #print(get(paste0('df4$', x)))
    }
    
    #print(BadCol)  
    #remove predictors over cut threshold 
    bad <- names(data) %in% c(BadCol, remove)
    good <- data[!bad]
    #remove predictors with NA
    #removeNA <- names(good) %in% noNA
    #names <- good[removeNA]
    names <- good
    names <- names[,9:ncol(names)]

    # #count the TRUE
    # df5 <- apply(df4, 1, function(x) length(which(x=="TRUE")))
    # df6 <- data.frame(t(sapply(df5,c)))
    # #remove any which appear less than twice (remember all appear once for correlating to themselves)
    # tempFunc <- function(VEC) max(VEC, na.rm  = TRUE) <= 1 | any(is.na(VEC))
    # BadCol <- apply(X = df6, MARGIN = 2,  tempFunc)
    # names <- df6[,!BadCol]
    
    #subset the data
    if(nchar(preserve[1])>0){
    data <- data[, (colnames(data) %in% c('sppres','x','y',colnames(names), preserve))]
        }else{
    data <- data[, (colnames(data) %in% c('sppres','x','y',colnames(names)))]
        }
    #data = na.omit(data)
    #uncomment to check for and remove NAs
    data = data[ , colSums(is.na(data)) == 0]
    output <- list(data,df1)
} 

boxFunction <- function(cutData){
  name <- names(cutData[,4:ncol(cutData)])
  for(n in name) {
    x1 <- unlist(cutData[n])
    x2 <- unlist(cutData['sppres'])
    boxplot(x1~x2,main = n, xlab = "Absence/Presence", ylab = n)
  }
}

#This function is designed to make the raster stack and subset by the selected columns
rasterStack <- function(data, rasterList, column_names) {   
    newList <- paste0("/vsicurl/https://storage.googleapis.com/predictors_public/",rasterList)
    newr <- stack(newList)
    names(newr) <- column_names
    cut <- subset(newr, colnames(data[,4:ncol(data)]))
    cut
    }

#Generalized linear model
library(PresenceAbsence)
library(DAAG)
library(grid)
library(gridExtra)
formatColumns <- function(dat) {
    resp <- colnames(dat[1]) # assign resp column name
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    mod.formula <- as.formula(paste(as.factor(resp), "~", paste(pred, collapse = "+"))) # formula
}
glmFunction <- function(cutData, rasters) {
  prestep <- glm(formatColumns(cutData), family = binomial, data = cutData)
  model <-step(prestep, trace = F)
  prediction <- predict(model, type = "response")
  mod1 <- "model"
  data <-cbind(mod1, cutData[1], prediction)
  thresh <- optimal.thresholds(data, opt.methods = c("Default"))
  ##table
  # plot.new()
  table <- summary(model)
  t <- as.data.frame(table$coefficients) %>% 
    mutate_if(is.numeric, ~round(., 5))
  t <- t[order( t[,ncol(t)] ),]
  t <- head(t,20)
  # grid.draw(tableGrob(t))
  grid.newpage()
  vp <- viewport(x = 0.4, y = 0.35, width = 1, height = 5) 
  grid.rect(vp = vp)
  tg <- tableGrob(t, vp = vp)
  # tg$widths[-1] <- rep(unit(1/2,"null"), 2)
  # tg$heights <- rep(unit(1/nrow(tg), "null"), nrow(tg))
  grid.draw(tg)
  
  #start accuracy assessment

  nfolds <- nrow(cutData)
  jack <- CVbinary(model, nfolds = nfolds, rand = c(1:nfolds), print.details = F)
  jack <- jack$cvhat
  data2 <- cbind(data, jack)
  auc.roc.plot(data2, color = T)
  
  accuracy <- presence.absence.accuracy(data2, threshold = thresh$prediction, st.dev = F)
  tss <- accuracy$sensitivity+accuracy$specificity-1
  
  accuracy <- cbind(accuracy[1:7], tss) %>% mutate_if(is.numeric, ~round(., 5))
  plot.new()
  grid.draw(tableGrob(accuracy))
  #end acc assess
  prediction2 <- predict(rasters, model, type = "response", fun = predict, index = 2)
  reclass <- reclassify(prediction2, (c(0, thresh$prediction, 0, thresh$prediction, 1, 1)))
  reclass
}

#Generalized additive model
library(gam)
formatGam <- function(dat) {
    resp <- colnames(dat[1]) # assign resp column name
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    form <- as.formula(paste(as.factor(resp), "~", paste(paste("lo(",pred,",5)"), collapse = "+"))) # formula
}

formatScope <- function(dat) {
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    form2 <- as.list(paste("~1 +",pred,"+ lo(",pred,", 3) + lo(",pred,",5)")) #
}
gamFunction <- function(cutData, rasters) {
  prestep <- gam(formatGam(cutData), family = binomial, data = cutData)
  out<- formatScope(cutData) #new
  out<- lapply(out,as.formula) ##
  names(out) <- colnames(cutData[c(4:ncol(cutData))]) # new
  model <-step.Gam(prestep, scope=out) #new
  #table
  table <- summary(model)
  t <- as.data.frame(table$anova) %>%
    mutate_if(is.numeric, ~round(., 5)) %>% 
    na.omit()
  t <- t[order( t[,1] ),]
  grid.newpage()
  vp <- viewport(x = 0.4, y = 0.18, width = 1, height = 5) 
  grid.rect(vp = vp)
  tg <- tableGrob(t, vp = vp)
  grid.draw(tg)
  #
  prediction <- predict(model, type = "response")
  mod1 <- "model"
  data <-cbind(mod1, cutData[1], prediction)
  thresh <- optimal.thresholds(data, opt.methods = c("Default"))
  #
  nfolds <- nrow(cutData)
  jack <- CVbinary(model, nfolds = 3, print.details = F)
  jack <- jack$cvhat
  data2 <- cbind(data, jack)
  auc.roc.plot(data2, color = T)
  #
  #start acc ass
  accuracy <- presence.absence.accuracy(data2, threshold = thresh$prediction, st.dev = F) #data
  tss <- accuracy$sensitivity+accuracy$specificity -1
  accuracy <- cbind(accuracy[1:7], tss) %>% mutate_if(is.numeric, ~round(., 5))
  plot.new()
  grid.draw(tableGrob(accuracy))
  
  #
  prediction2 <- predict(rasters, model, type = "response", fun = predict, index = 2)
  reclass <- reclassify(prediction2, (c(0, thresh$prediction, 0, thresh$prediction, 1, 1)))
  reclass 
}


#MaxEnt
library(dismo)
library(rJava)
library(data.table)#as.data.table(
maxFunction <- function(cutData, rasters) {   
  xy <- cutData[,2:3]
  fold <- kfold(xy, k = 5)
  validation <- xy[fold == 1, ]
  training <- xy[fold!= 1, ]
  
  model <- maxent(rasters, training)
  randomp <- randomPoints(rasters, 1000)
  plot(model,cex=.5)
  
  modelxy <- data.frame(randomPoints(rasters, dim(xy)[1]*4))
  modelpts <- data.frame(raster::extract(rasters, modelxy))
  modelxy$sppres <- as.integer(0)
  
  prespts <- data.frame(raster::extract(rasters, xy))
  prespred <- predict(model, prespts)
  eval <- evaluate(model, p = validation, a = modelxy[,-3], x = rasters)
  truepres <- cutData %>% filter(sppres == 1)
  
  modl <- "model"
  
  modelprediction1 <- cbind(as.data.table(modl), as.data.table(truepres[1]), as.data.table(prespred))
  prespred <- predict(model, modelpts)
  #pad shortest vector with NA's to have same length as longest vector
  length(prespred) <- length(modelxy[3])
  #
  modelprediction2 <- cbind(modl, modelxy[3], prespred)
  data <- data.frame(rbind(modelprediction1, modelprediction2))
  threshold <- threshold(eval)
  #optimal.thresholds(data, opt.methods = 1:6)
  #threshold[c(1:2)]
  #optimal.thresholds(data, opt.method = c(4,3))
  #mod1.cfmat <- table(data[[2]], factor(as.numeric(data$prespred >= threshold$spec_sens)))
  #acc ass
  auc.roc.plot(data, color = T)
  accuracy <- presence.absence.accuracy(data, threshold = threshold$spec_sens, st.dev = F)
  tss <- accuracy$sensitivity + accuracy$specificity - 1
  accuracy <- cbind(accuracy[1:7], tss) %>% mutate_if(is.numeric, ~round(., 5))
  plot.new()
  grid.draw(tableGrob(accuracy))
  #
  modelprob = predict(model, rasters)
  modelclas = reclassify(modelprob, c(0,threshold[[2]],0,threshold[[2]],1,1))
  modelclas 
  
} 

# -
#Random Forest model
library(randomForest)
formatRaf <- function(dat) {
    resp <- colnames(dat[1]) # assign resp column name
    resp <- paste("as.factor(", colnames(dat[1]), ")", sep = "")
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    mod.formula <- as.formula(paste(resp, "~", paste(pred, collapse = "+"))) # formula
}
rafFunction <- function(cutData, rasters) {   
    model <- randomForest(formatRaf(cutData), importance = T, keep.forest = T, data = cutData)
    varImpPlot(model,type=2)
    mod1.pred <- predict(model, type = "prob")[,2]
    modl <- "model"
    data <- cbind(modl, cutData[1], mod1.pred)
    auc.roc.plot(data, color = T)
    #ac as
    threshold <- optimal.thresholds(data, opt.methods = c("ReqSens"), req.sens = 0.95)
    accuracy <- presence.absence.accuracy(data, threshold = threshold$mod1.pred, st.dev = F)
    tss <- accuracy$sensitivity + accuracy$specificity -1
    accuracy <- cbind(accuracy[1:7], tss) %>% mutate_if(is.numeric, ~round(., 5))
    plot.new()
    grid.draw(tableGrob(accuracy)) 
    #
    # oob.acc <- presence.absence.accuracy(data, st.dev = F)
    # tss <-oob.acc$sensitivity + oob.acc$specificity - 1
    # oob.acc <- cbind(oob.acc[1:7], tss)
    # oob.acc$model <- "mod1.predoob"
    # RF.acc <- rbind(oob.acc[c(1, 4:5, 7:8)],accuracy[c(1, 4:5, 7:8)])
    modelprob <- predict(rasters, model, type="response", fun = predict, index = 2)
    modelclas <- reclassify(modelprob, c(0,threshold[[2]],0,threshold[[2]],1,1))
    modelclas
    }


#Boosted Regression Tree model
library(gbm)
brtFunction <- function(data, rasters) {   
    dat1 <- na.omit(data)
    n.col <- ncol(dat1)
    pred <- 4:n.col
    resp <- paste("as.factor(", colnames(dat1[1]), ")", sep = "")
    model <- gbm.step(data = dat1, gbm.x = pred, gbm.y = 1, family = "bernoulli", tree.complexity = 3, learning.rate = 1e-04, bag.fraction = .75, n.folds = 10, n.trees = 50, plot.main = TRUE, keep.fold.fit = T)
    #mod2.BRT <-gbm.step(data = dat1, gbm.x = pred, gbm.y = 1, family = "bernoulli", tree.complexity = 3, learning.rate = .1, bag.fraction = .75, n.folds = 10, plot.main = TRUE, keep.fold.fit = TRUE)    
    modl <- "mod2.BRT"
    dat2 <- cbind(modl, dat1[1], model$fitted, model$fold.fit)
    names(dat2)[3:4] <-c("pred", "cvpred")
    dat2$cvpred <- exp(dat2$cvpred)/(1+ exp(dat2$cvpred))
    threshold <- optimal.thresholds(dat2, opt.methods = c("ObsPrev"))
    #mod1.cfmatR <- table(dat2[[2]], factor(as.numeric(dat2$pred >= threshold$pred)))
    gbm.plot(model)

    #mod1.cfmatX <- table(dat2[[2]], factor(as.numeric(dat2$cvpred >= threshold$cvpred)))
    #acc as
    auc.roc.plot(dat2, color = T)
    accuracy <- presence.absence.accuracy(dat2, threshold = threshold$pred, st.dev = F)
    tss <- accuracy$sensitivity + accuracy$specificity - 1
    accuracy.brt <- cbind(accuracy[1:7], tss)  %>% mutate_if(is.numeric, ~round(., 5))
    plot.new()
    grid.draw(tableGrob(accuracy.brt))
    #
    modelprob = predict(rasters, model, n.trees = model$gbm.call$best.trees, type = "response")
    modelclas <- reclassify(modelprob, c(0,threshold[[2]],0,threshold[[2]],1,1))
    
    modelclas
}
