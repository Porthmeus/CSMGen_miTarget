# Porthmeus
# 16.12.21

dbg <- FALSE
if(!dbg){
    log <- file(snakemake@log[[1]], open = "wt")
    sink(log)
    sink(log, type = "message")
}

# set max-ppsize  - whatever that is...
options(expressions = 50000)
# load libraries
require(caret)
require(data.table)


# load data
if(dbg){
    mat.rxns <- read.csv("results/ModelAnalysis/rxnIdMat.GL75.MatjesAbsorption.colormore3D.csv", row.names = 1)
    mat.center <- read.csv("results/FVA/centerFVA.GL75.MatjesAbsorption.colormore3D.csv", row.names = 1)
    mat.range <- read.csv("results/FVA/rangeFVA.GL75.MatjesAbsorption.colormore3D.csv", row.names = 1)
    # get the meta data and merge them
    meta <- fread("./resources/META_emed_future.csv")
    clinic <- fread("./resources/ClinicalData.csv")
} else {
    mat.rxns <- read.csv(snakemake@input[["rxns"]], row.names = 1)
    mat.center <- read.csv(snakemake@input[["fva_center"]], row.names = 1)
    mat.range <- read.csv(snakemake@input[["fva_range"]], row.names = 1)
    meta <- fread(snakemake@input[["meta"]])
    clinic <- fread(snakemake@input[["clinic"]])
}

mat.rxns2 <- matrix(0, ncol = ncol(mat.rxns), nrow = nrow(mat.rxns), dimnames = list(rownames(mat.rxns), colnames(mat.rxns)))
mat.rxns2[mat.rxns == "True"] <- 1
mat.rxns <- mat.rxns2
meta[,PatientID_Time := paste(PatientID, Time_seq, sep = "_")]
meta <- merge(meta, clinic, by = "PatientID_Time")

# merge data
mat.all <- merge(t(mat.range), t(mat.center), by = 0, suffix = c(".FVA_range",".FVA_center"))
mat.all <- merge(mat.all, t(mat.rxns2), by.x = "Row.names",by.y = 0, suffix = c("",".rxns"))
mat.all <- merge(meta[,.(SeqID, HB_Mayo_impu)], mat.all, by.x = "SeqID", by.y = "Row.names")

# drop the ID column
rn <- mat.all[,SeqID] 
mat.all <- as.data.frame(mat.all[,2:ncol(mat.all)])
rownames(mat.all) <- rn

## preprocessing
# remove columns without variance
total <- ncol(mat.all)
vars <- nearZeroVar(mat.all)
#vars <- apply(mat.all, 2, var)
print(paste0("Will remove ",
             length(vars),
             " (",
             round(length(vars)/total*100,digits=2),
             "%) columns because of near zero variance!"))
mat.all <- mat.all[,-vars]

# remove columns with NAs
noNAs <- !apply(mat.all, 1, function(x) any(is.na(x)))
mat.all <- mat.all[noNAs,]
print(paste0("Removed ", sum(!noNAs), " samples due to NAs"))

# remove linear combinations
#comboInfo <- findLinearCombos(mat.all)
#mat.all.ori <- mat.all
#mat.all <- mat.all[,-comboInfo[["remove"]]]



# add the ML algorithms to test
mls <- c("ranger", # random forest
         "knn" # K-nearest neighbor
         #"bayesglm" # bayesian generalized linear model
         #"svmPoly", # support vector machine with polynomial kernel
         #"enet" # elestic net
)

# define how often to run the training and test it against a test set
if(dbg){
    # create a smaller subset to test the code
    cols <- c(colnames(mat.all)[1],sample(colnames(mat.all)[-1], 40))
    mat.all <- mat.all[,cols]
    # test mls
    mls <- c("knn","enet")
    n <- 5
} else {
    n <- snakemake@params[["n_repeats"]]
}

# define variables to store the data in
stats <- list(
              "model" = rep("", n*length(mls)),
              "R2" = rep(0,n*length(mls)),
              "RMSE" = rep(0, n*length(mls)),
              "bias" = rep(0, n*length(mls))
              )
idx <- 0
for( i in 1:n){
    for( j in 1:length(mls)){
        ml <- mls[j]
        idx <- idx+1

        # partition data to test and to training data
        train = unlist(createDataPartition(mat.all$HB_Mayo_impu, p = 0.75, list = TRUE))
        test <- rownames(mat.all)[-train]
        train <- rownames(mat.all)[train]

        # consider the independence of sampling due to time series in the resampling approach
        train.patients <- meta[SeqID %in% train, PatientID.x]
        train.patients.folds <- groupKFold(train.patients, k = 5)

        # scale data
        outcome.train <- mat.all[train,1]
        outcome.test <- mat.all[test,1]
        scaler <- preProcess(mat.all[train,-1], method = c("center","scale"))
        training <- predict(scaler, mat.all[train,-1])
        
        # build a ml model
        fit.rf <- train(#data = training,
                        #data = mat.all[train,],
                        #HB_Mayo_impu~.,
                        x = training,
                        y = outcome.train,
                        method = ml,
                        trControl = trainControl(method = "LGOCV", number = 5, search = "random", index = train.patients.folds),
                        tuneLength = 10)

        # test the model and store the metrics in the list
        testing <- predict(scaler, mat.all[test,])
        test.predict <- predict(fit.rf, testing)
        stats[["model"]][idx] <- ml
        stats[["R2"]][idx] <-cor(test.predict, outcome.test)
        stats[["RMSE"]][idx]<-sd(test.predict- outcome.test)
        stats[["bias"]][idx]<-mean(test.predict) - mean(outcome.test) 

        # report progress
        timestamp()
        print(paste0(idx, " of ", length(mls)*n, " (", round(idx/(n*length(mls))*100), "%) models build"))
    }
}

stats <-data.table(as.data.frame(stats))

if(dbg){
    write.csv(stats, file = "temp/MLstats.test.csv", row.names = FALSE)
} else {
    write.csv(stats, file = snakemake@output[["stats"]], row.names = FALSE)
}
