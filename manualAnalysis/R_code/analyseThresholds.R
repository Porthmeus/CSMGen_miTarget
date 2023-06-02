# Porthmeus
# 21.12.21

## @knitr loadLibraries
require(data.table)
require(ggplot2)
require(RColorBrewer)
require(DT)


## @knitr loadData
# get all the files
dr <- "../results/testCoreThresholds"
fls <- list.files(dr)
fls <- file.path(dr,fls)

# go through the files and create a large data frame with the meta information
rm("tbl.MLstats.all")
for(fl in fls){
    info <- strsplit(basename(fl), split = ".", fixed = TRUE)[[1]]
    tbl.MLstats <- read.csv(fl)
    tbl.MLstats <- cbind(tbl.MLstats, Threshold = info[2])
    tbl.MLstats <- cbind(tbl.MLstats, Diet = info[3])
    tbl.MLstats <- cbind(tbl.MLstats, MetModel = info[4])
    if(exists("tbl.MLstats.all")){
        tbl.MLstats.all <- rbind(tbl.MLstats.all, tbl.MLstats)
    } else {
        tbl.MLstats.all <- tbl.MLstats
    }
}

tbl.MLstats.all <- data.table(tbl.MLstats.all)

## @knitr R2plot
pp_R2 <- ggplot(tbl.MLstats.all, aes(x = Threshold, y = R2, fill = Threshold)) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(model~MetModel)+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))+
    scale_fill_brewer(palette = "Dark2")
pp_R2

## @knitr RMSEplot
pp_RMSE <- ggplot(tbl.MLstats.all, aes(x = Threshold, y = RMSE, fill = Threshold)) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(model~MetModel)+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))+
    scale_fill_brewer(palette = "Dark2")
pp_RMSE

## @knitr biasPlot
pp_bias <- ggplot(tbl.MLstats.all, aes(x = Threshold, y = bias, fill = Threshold)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(model~MetModel)+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))+
    scale_fill_brewer(palette = "Dark2")
pp_bias

## @knitr scorePlot
tbl.MLstats.all[,score := R2/((RMSE+sqrt(bias^2))*0.5)]

pp_score <- ggplot(tbl.MLstats.all, aes(x = Threshold, y = score, fill = Threshold)) +
    geom_boxplot() +
    theme_bw() +
    facet_wrap(~MetModel)+
    theme(axis.text.x = element_text(angle = 45, hjust=1, vjust = 1))+
    scale_fill_brewer(palette = "Dark2")
pp_score

## @knitr summaryTable
tbl.MLstats.model <- tbl.MLstats.all[, .(R2 = mean(R2, na.rm = TRUE),
                                         RMSE = mean(RMSE, na.rm = TRUE),
                                         bias = mean(bias, na.rm = TRUE),
                                         score = mean(score, na.rm = TRUE),
                                         R2_sd = sd(R2, na.rm = TRUE),
                                         RMSE_sd = sd(RMSE, na.rm = TRUE),
                                         bias_sd = sd(bias, na.rm = TRUE),
                                         score_sd = sd(score, na.rm = TRUE),
                                         bias_sem = sd(bias, na.rm = TRUE)/sqrt(.N),
                                         R2_sem = sd(R2, na.rm = TRUE)/sqrt(.N),
                                         RMSE_sem = sd(RMSE, na.rm = TRUE)/sqrt(.N),
                                         score_sem = sd(score, na.rm = TRUE)/sqrt(.N)),
                                by = .(MetModel, Diet, model, Threshold)]
tbl.MLstats.model[score == max(score),]

## @knitr bestByLMcmplx
mod <- lm(data = tbl.MLstats.all, score ~ 0+Threshold*MetModel)
mod.anova <- anova(mod)
datatable(as.data.frame(mod.anova))
mod.sum <- as.data.frame(coef(summary(mod)))
datatable(mod.sum)

## @knitr bestByLMsimple
mod <- lm(data = tbl.MLstats.all, score ~ 0+Threshold:MetModel)
mod.sum <- as.data.frame(coef(summary(mod)))
datatable(mod.sum)
