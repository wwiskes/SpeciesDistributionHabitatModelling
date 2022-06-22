# ---
# This script is a companion to the main script titled "SDHMmain"
# This script stores functions called by the main script and should not be altered.
# Written by William Wiskes
# Last update 6/21/2021
# ---

#this function will pull in data from a postgres database
library(RPostgreSQL)
queryPostgres <- function(species) {
    source("env.R") #set all the variables needed from a postgres database in an external env file
    con <- dbConnect(RPostgres::Postgres(), dbname = db, host=host_db, port=db_port, user=db_user, password=db_password)  
    query <- dbSendQuery(con, paste0("SELECT * FROM ",species))
    occs <- dbFetch(query)
    occs %>% drop_na()
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
        d <- cbind(d, extract(c, pointData[,c("x","y")]))
        i = i + 1
    }
d
}

#To generate pseudo absences from the extent of the study area
library(sfheaders)
pseudoFunction <- function(pointSF, blob, temp) {
    #intersect points with state, this is redundant with this data set, but may be important for other workflows
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

    #Make fnet into dataframe

    fnetSF <- st_crop(fnetSF, blob)
    fnetDF <- as.data.frame(fnetSF)
    fnetDF <- fnetDF[,1:3]

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

#This function removes columns from the data without statistical relevance according to spearmans 
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
    #return a dataframe which all that have a correlation above the abs of cut as TRUE
    df4 <-as.data.frame(subset(df3 > cut | df3 < -cut))
    BadCol <- c()
    for (x in colnames(df4)) {
      var <- df4 %>%
        select(x) %>%
        filter(get(x) != FALSE)
      if (row.names(var) != colnames(var)) {
        #print(row.names(var)[2:length(row.names(var))])
        BadCol <- append(BadCol,row.names(var)[2:length(row.names(var))])
      }
      
      #print(get(paste0('df4$', x)))
    }
    #print(BadCol)  
    bad <- names(data) %in% c(BadCol, remove)
    names <- data[!bad]
    
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
    output <- list(data,df1)
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
formatColumns <- function(dat) {
    resp <- colnames(dat[1]) # assign resp column name
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    mod.formula <- as.formula(paste(as.factor(resp), "~", paste(pred, collapse = "+"))) # formula
}
glmFunction <- function(data, cut) {
    mod1.LR <- glm(formatColumns(data), family = binomial, data = data)
    mod2.LR <-step(mod1.LR, trace = F)
    mod2.pred <- predict(mod2.LR, type = "response")
    mod1 <- "mod2.LR"
    dat2 <-cbind(mod1, data[1], mod2.pred)
    mod.cut.GLM <- optimal.thresholds(dat2, opt.methods = c("Default"))
    
    # newList <- gsub("gs://", "/vsicurl/https://storage.googleapis.com/",rasterList)
    # # newKeep <- gsub("gs://", "/vsicurl/https://storage.googleapis.com/",rasterKeep)
    # # appended <- append(newList, newKeep)
    # newr <- stack(newList)
    # names(newr) <- column_names
    # cut <- subset(newr, colnames(data[,4:ncol(data)]))

    modFprob.LR <- predict(cut, mod2.LR, type = "response", fun = predict, index = 2)
    modFclas.LR <- reclassify(modFprob.LR, (c(0, mod.cut.GLM$mod2.pred, 0, mod.cut.GLM$mod2.pred, 1, 1)))
    modFclas.LR
}

#Generalized additive model
library(gam)
formatGam <- function(dat) {
    resp <- colnames(dat[1]) # assign resp column name
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    mod.form <- as.formula(paste(as.factor(resp), "~", paste(paste("lo(",pred,",5)"), collapse = "+"))) # formula
}
# formatScope <- function(dat) {
#     resp <- colnames(dat[1]) # assign resp column name
#     pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
#     mod.form2 <- as.formula(paste(pred,"= ~1 +",pred,"+ lo(",pred,", 3) + lo(",pred,",5)")) # remove as.formula
# }
formatScope <- function(dat) {
    pred <- colnames(dat[c(4:ncol(dat))]) # assign preds column names
    mod.form2 <- as.list(paste("~1 +",pred,"+ lo(",pred,", 3) + lo(",pred,",5)")) #
}
gamFunction <- function(data, cut) {
    mod1.LR <- gam(formatGam(data), family = binomial, data = data)
    out<- formatScope(data) #new
    out<- lapply(out,as.formula) ##
    names(out) <- colnames(data[c(4:ncol(data))]) # new
    mod2.LR <-step.Gam(mod1.LR, scope=out) #new
    #mod2.LR <-step.Gam(mod1.LR, scope=list(formatScope(data))) # old
    mod2.pred <- predict(mod2.LR, type = "response")
    mod1 <- "mod2.LR"
    dat2 <-cbind(mod1, data[1], mod2.pred)
    mod.cut.GLM <- optimal.thresholds(dat2, opt.methods = c("Default"))
    
    # newList <- gsub("gs://", "/vsicurl/https://storage.googleapis.com/",rasterList)
    # # newKeep <- gsub("gs://", "/vsicurl/https://storage.googleapis.com/",rasterKeep)
    # # appended <- append(newList, newKeep)
    # newr <- stack(newList)
    # names(newr) <- column_names
    # cut <- subset(newr, colnames(data[,4:ncol(data)]))

    modFprob.LR <- predict(cut, mod2.LR, type = "response", fun = predict, index = 2)
    modFclas.GAM <- reclassify(modFprob.LR, (c(0, mod.cut.GLM$mod2.pred, 0, mod.cut.GLM$mod2.pred, 1, 1)))
    modFclas.GAM 
}


#MaxEnt model
library(dismo)
library(rJava)

maxFunction <- function(data, cut) {   
    xy <- data[,2:3]
    fold <- kfold(xy, k = 5)
    pres.tst <- xy[fold == 1, ]
    pres.tr <- xy[fold!= 1, ]
    
    mod1.MAX <- maxent(cut, pres.tr)
    mod2.MAX <- mod1.MAX
    mod2.bak <- randomPoints(cut, 1000)

    mod2.val <- evaluate(mod2.MAX, p = pres.tst, a = mod2.bak, x = cut)
    pts.tst <- data.frame(extract(cut, pres.tst))
    pts.bak <- data.frame(extract(cut, mod2.bak))
    mod2.val <- evaluate(mod2.MAX, p = pts.tst, a = pts.bak)

    bak.xy <- data.frame(randomPoints(cut, dim(xy)[1]*4))
    bak.pts <- data.frame(extract(cut, bak.xy))
    bak.pred <- predict(mod1.MAX, bak.pts)
    bak.xy$sppres <- as.integer(0)

    pres.pts <- data.frame(extract(cut, xy))
    pres.pred <- predict(mod1.MAX, pres.pts)
    mod1.val <- evaluate(mod1.MAX, p = pres.tst, a = bak.xy[,-3], x = cut)
    ybcu.tr.p <- data %>% filter(sppres == 1)

    modl <- "mod1.MAX"
    mod1.pred <- pres.pred
    tmp.p <- cbind(modl, ybcu.tr.p[1], mod1.pred)
    mod1.pred <- bak.pred
    tmp.b <- cbind(modl, bak.xy[3], mod1.pred)
    dat2 <- data.frame(rbind(tmp.p, tmp.b))
    mod.cut.MAXENT <- threshold(mod1.val)
    optimal.thresholds(dat2, opt.methods = 1:6)
    mod.cut.MAXENT[c(1:2)]
    optimal.thresholds(dat2, opt.method = c(4,3))
    mod1.cfmat <- table(dat2[[2]], factor(as.numeric(dat2$mod1.pred >= mod.cut.MAXENT$spec_sens)))
    mod1.acc <- presence.absence.accuracy(dat2, threshold = mod.cut.MAXENT$spec_sens, st.dev = F)
    tss <- mod1.acc$sensitivity + mod1.acc$specificity - 1
    mod1.acc <- cbind(mod1.acc[1:7], tss)

    mod1.MAXprob = predict(mod1.MAX, cut)
    mod1.MAXclas = reclassify(mod1.MAXprob, c(0,mod.cut.MAXENT[[2]],0,mod.cut.MAXENT[[2]],1,1))
    mod1.MAXclas 
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
rafFunction <- function(data, cut) {   
    mod1.RF <- randomForest(formatRaf(data), importance = T, keep.forest = T, data = data)
    mod1.pred <- predict(mod1.RF, type = "prob")[,2]
    modl <- "mod1.RF"
    dat2 <- cbind(modl, data[1], mod1.pred)
    
    # newList <- gsub("gs://", "/vsicurl/https://storage.googleapis.com/",rasterList)
    # newr <- stack(newList)
    # names(newr) <- column_names
    # cut <- subset(newr, colnames(data[,4:ncol(data)]))
    
    mod.cut.RF <- optimal.thresholds(dat2, opt.methods = c("ReqSens"), req.sens = 0.95)
    mod1.acc <- presence.absence.accuracy(dat2, threshold = mod.cut.RF$mod1.pred, st.dev = F)
    tss <- mod1.acc$sensitivity + mod1.acc$specificity -1
    mod1.acc <- cbind(mod1.acc[1:7], tss)
    oob.acc <- presence.absence.accuracy(dat2, st.dev = F)
    tss <-oob.acc$sensitivity + oob.acc$specificity - 1
    oob.acc <- cbind(oob.acc[1:7], tss)
    oob.acc$model <- "mod1.predoob"
    RF.acc <- rbind(oob.acc[c(1, 4:5, 7:8)],mod1.acc[c(1, 4:5, 7:8)])
    mod1.RFprob <- predict(cut, mod1.RF, type="response", fun = predict, index = 2)
    mod1.RFclas <- reclassify(mod1.RFprob, c(0,mod.cut.RF[[2]],0,mod.cut.RF[[2]],1,1))
    mod1.RFclas
    }


#Boosted Regression Tree model
library(gbm)
brtFunction <- function(data, cut) {   
    dat1 <- na.omit(data)
    n.col <- ncol(dat1)
    pred <- 4:n.col
    resp <- paste("as.factor(", colnames(dat1[1]), ")", sep = "")
    mod1.BRT <- gbm.step(data = dat1, gbm.x = pred, gbm.y = 1, family = "bernoulli", tree.complexity = 3, learning.rate = 1e-04, bag.fraction = .75, n.folds = 10, n.trees = 50, plot.main = TRUE, keep.fold.fit = T)
    mod2.BRT <-gbm.step(data = dat1, gbm.x = pred, gbm.y = 1, family = "bernoulli", tree.complexity = 3, learning.rate = .1, bag.fraction = .75, n.folds = 10, plot.main = TRUE, keep.fold.fit = TRUE)    
    modl <- "mod2.BRT"
    dat2 <- cbind(modl, dat1[1], mod1.BRT$fitted, mod1.BRT$fold.fit)
    names(dat2)[3:4] <-c("pred", "cvpred")
    dat2$cvpred <- exp(dat2$cvpred)/(1+ exp(dat2$cvpred))
    mod.cut.BRT <- optimal.thresholds(dat2, opt.methods = c("ObsPrev"))
    mod1.cfmatR <- table(dat2[[2]], factor(as.numeric(dat2$pred >= mod.cut.BRT$pred)))

    mod1.cfmatX <- table(dat2[[2]], factor(as.numeric(dat2$cvpred >= mod.cut.BRT$cvpred)))
    mod1.acc <- presence.absence.accuracy(dat2, threshold = mod.cut.BRT$pred, st.dev = F)
    tss <- mod1.acc$sensitivity + mod1.acc$specificity - 1

    mod1.acc.brt <- cbind(mod1.acc[1:7], tss)
    mod1.BRTprob = predict(cut, mod1.BRT, n.trees = mod1.BRT$gbm.call$best.trees, type = "response")
    mod1.BRTclas <- reclassify(mod1.BRTprob, c(0,mod.cut.BRT[[2]],0,mod.cut.BRT[[2]],1,1))
    mod1.BRTclas
}
