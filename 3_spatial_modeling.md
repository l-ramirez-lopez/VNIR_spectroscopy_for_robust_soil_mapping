========================================================

Refers to section Spatial modelling of spectrally predicted soil
properties

------------------------------------------------------------------------

Install the required packages

    Sys.setenv(language = "EN")
    options(warn=-1)
    if (!require(resemble))
      install.packages("resemble") # version 1.2.2
    if (!require(prospectr))
      install.packages("prospectr") # version 0.1.3
    if (!require(clhs))
      install.packages("clhs") # version 0.5-6
    if (!require(matrixStats))
      install.packages("matrixStats") # version 0.51.0
    if (!require(doParallel))
      install.packages("doParallel") # version 1.0.10
    if (!require(ggplot2))
      install.packages("ggplot2") # version 2.2.1
    if (!require(caret))
      install.packages("caret") # version 6.0-73
    if (!require(georob))
      install.packages("georob") # version 0.3-3
    if (!require(rgdal))
      install.packages("rgdal") # version 1.2-7
    if (!require(sp))
      install.packages("sp") # version 1.2-4
    if (!require(raster))
      install.packages("raster") # version 2.5-8
    if (!require(gridExtra))
      install.packages("gridExtra") # version 2.2.1
    if (!require(RColorBrewer))
      install.packages("RColorBrewer") # version 1.1-2
    if (!require(grid))
      install.packages("grid") # version 3.3.2
    if (!require(cowplot))
      install.packages("cowplot") # version 1.1-2
    if (!require(extrafont))
      install.packages("extrafont") # version 3.3.2
    # for the first time only, run:
    ## font_import()
    ## loadfonts(device = "win")
    options(warn=0)

Define the working directory

    workingd <- "your/working/directory/"
    setwd(workingd)

Organize the data
=================

    data.spec <- read.table("spec_aideddb.txt", header = TRUE, check.names = FALSE, sep ="\t")

Load the variance of the spectroscopic model error

    load('var_spec_error.Rdata')

Log-transformation of the target variables

    data.spec$alr_Silt_lab <- log(data.spec[,"Silt"]/data.spec[,"Sand"])
    data.spec$alr_Clay_lab <- log(data.spec[,"Clay"]/data.spec[,"Sand"])

    data.spec$alr_Silt_spec <- log(data.spec[,"specSilt_alr"]/data.spec[,"specSand_alr"])
    data.spec$alr_Clay_spec <- log(data.spec[,"specClay_alr"]/data.spec[,"specSand_alr"])

Make some summary statistics of the target variables

    summary_stats_a <- rbind(t(colQuantiles(as.matrix(data.spec[data.spec$layer == "A", 
                                                          c("Sand", "Silt", "Clay", "Ca")]))),
                             colMeans(data.spec[data.spec$layer == "A", 
                                                          c("Sand", "Silt", "Clay", "Ca")]),
                             colSds(as.matrix(data.spec[data.spec$layer == "A", 
                                                          c("Sand", "Silt", "Clay", "Ca")])))

    summary_stats_a <- round(summary_stats_a, 2)
    rownames(summary_stats_a) <- c("Minimum", "1st Quantile", "Median", 
                                   "3rd quantile", "Maximum", "Mean", "Standard deviation")

    summary_stats_b <- rbind(t(colQuantiles(as.matrix(data.spec[data.spec$layer == "B", 
                                                          c("Sand", "Silt", "Clay", "Ca")]))),
                             colMeans(data.spec[data.spec$layer == "B", 
                                                          c("Sand", "Silt", "Clay", "Ca")]),
                             colSds(as.matrix(data.spec[data.spec$layer == "B", 
                                                          c("Sand", "Silt", "Clay", "Ca")])))

    summary_stats_b <- round(summary_stats_b, 2)
    rownames(summary_stats_b) <- c("Minimum", "1st Quantile", "Median", 
                                   "3rd quantile", "Maximum", "Mean", "Standard deviation")


    summary_stats <- rbind(data.frame(Property = colnames(summary_stats_b), 
                                      t(summary_stats_a), layer = "A"), 
                           data.frame(Property = colnames(summary_stats_b), 
                                      t(summary_stats_b), layer = "B"))

    summary_stats
    write.table(summary_stats, file = "summary_stats_BB.txt", sep = "\t", row.names = F)

fitting isotropic model
-----------------------

---- Models layer A ----
========================

alr\_Silt\_lab

    vario_alr_Silt_lab_a <- sample.variogram(data.spec[data.spec$layer == "A" & data.spec$set != "validation","alr_Silt_lab"], 
                                             locations = data.spec[data.spec$layer == "A" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                             lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Silt_lab_a[,c("lag.dist")], vario_alr_Silt_lab_a[,c("gamma")])

    alr_Silt_lab_a_psi2 <- georob(alr_Silt_lab ~ 1, data = data.spec[data.spec$layer == "A" & data.spec$set != "validation",], 
                                   locations = ~ POINT_X + POINT_Y, 
                                   variogram.model = "RMexp",
                                   param = c(variance = 0.1, nugget = 0.3, scale = 1500),
                                   fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                   tuning.psi = 2)


    summary(alr_Silt_lab_a_psi2)

alr\_Silt\_spec

    vario_alr_Silt_spec_a <- sample.variogram(data.spec[data.spec$layer == "A" & data.spec$set != "validation","alr_Silt_spec"], 
                                              locations = data.spec[data.spec$layer == "A"  & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                              lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Silt_spec_a[,c("lag.dist")], vario_alr_Silt_spec_a[,c("gamma")])

    alr_Silt_spec_a_psi2 <- georob(alr_Silt_spec ~ 1, data = data.spec[data.spec$layer == "A" & data.spec$set != "validation",], 
                                    locations = ~ POINT_X + POINT_Y, 
                                    variogram.model = "RMexp",
                                    param = c(variance = 1, nugget = 0.1, scale = 1500),
                                    fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                    tuning.psi = 2)


    summary(alr_Silt_spec_a_psi2)

alr\_Clay\_lab

    vario_alr_Clay_lab_a <- sample.variogram(data.spec[data.spec$layer == "A" & data.spec$set != "validation","alr_Clay_lab"], 
                                             locations = data.spec[data.spec$layer == "A"  & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                             lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Clay_lab_a[,c("lag.dist")], vario_alr_Clay_lab_a[,c("gamma")])

    alr_Clay_lab_a_psi2 <- georob(alr_Clay_lab ~ 1, data = data.spec[data.spec$layer == "A" & data.spec$set != "validation",], 
                                   locations = ~ POINT_X + POINT_Y, 
                                   variogram.model = "RMexp",
                                   param = c(variance = 1, nugget = 0.3, scale = 1500),
                                   fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                   tuning.psi = 2)

    summary(alr_Clay_lab_a_psi2)

alr\_Clay\_spec

    vario_alr_Clay_spec_a <- sample.variogram(data.spec[data.spec$layer == "A" & data.spec$set != "validation","alr_Clay_spec"], 
                                              locations = data.spec[data.spec$layer == "A" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                              lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Clay_spec_a[,c("lag.dist")], vario_alr_Clay_spec_a[,c("gamma")])

    alr_Clay_spec_a_psi2 <- georob(alr_Clay_spec ~ 1, data = data.spec[data.spec$layer == "A" & data.spec$set != "validation",], 
                                    locations = ~ POINT_X + POINT_Y, 
                                    variogram.model = "RMexp",
                                    param = c(variance = 1, nugget = 0.1, scale = 1500),
                                    fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                    tuning.psi = 2)

    summary(alr_Clay_spec_a_psi2)

Ca lab

    vario_Ca_lab_a <- sample.variogram(data.spec[data.spec$layer == "A" & data.spec$set != "validation","Ca"], 
                                       locations = data.spec[data.spec$layer == "A" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                       lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_Ca_lab_a[,c("lag.dist")], vario_Ca_lab_a[,c("gamma")])

    Ca_lab_a_psi2 <- georob(Ca ~ 1, data = data.spec[data.spec$layer == "A" & data.spec$set != "validation",], 
                             locations = ~ POINT_X + POINT_Y, 
                             variogram.model = "RMexp",
                             param = c(variance = 200, nugget = 10, scale = 500),
                             fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                             tuning.psi = 2)

    summary(Ca_lab_a_psi2)

specCa

    vario_specCa_a <- sample.variogram(data.spec[data.spec$layer == "A" & data.spec$set != "validation","specCa"], 
                                       locations = data.spec[data.spec$layer == "A" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                       lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_specCa_a[,c("lag.dist")], vario_specCa_a[,c("gamma")])

    specCa_a_psi2 <- georob(specCa ~ 1, data = data.spec[data.spec$layer == "A" & data.spec$set != "validation",], 
                             locations = ~ POINT_X + POINT_Y, 
                             variogram.model = "RMexp",
                             param = c(variance = 220, nugget = 1, scale = 1000),
                             fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                             tuning.psi = 2)


    summary(specCa_a_psi2)

---- Models layer B ----
========================

alr\_Silt\_lab

    vario_alr_Silt_lab_b <- sample.variogram(data.spec[data.spec$layer == "B" & data.spec$set != "validation","alr_Silt_lab"], 
                                             locations = data.spec[data.spec$layer == "B" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                             lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Silt_lab_b[,c("lag.dist")], vario_alr_Silt_lab_b[,c("gamma")])

    alr_Silt_lab_b_psi2<- georob(alr_Silt_lab ~ 1, data = data.spec[data.spec$layer == "B" & data.spec$set != "validation",], 
                                   locations = ~ POINT_X + POINT_Y, 
                                   variogram.model = "RMexp",
                                   param = c(variance = 0.1, nugget = 0.3, scale = 1500),
                                   fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                   tuning.psi = 2)

    summary(alr_Silt_lab_b_psi2)

alr\_Silt\_spec

    vario_alr_Silt_spec_b <- sample.variogram(data.spec[data.spec$layer == "B" & data.spec$set != "validation","alr_Silt_spec"], 
                                              locations = data.spec[data.spec$layer == "B"  & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                              lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Silt_spec_b[,c("lag.dist")], vario_alr_Silt_spec_b[,c("gamma")])

    alr_Silt_spec_b_psi2 <- georob(alr_Silt_spec ~ 1, data = data.spec[data.spec$layer == "B" & data.spec$set != "validation",], 
                                    locations = ~ POINT_X + POINT_Y, 
                                    variogram.model = "RMexp",
                                    param = c(variance = 1, nugget = 0.1, scale = 1500),
                                    fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                    tuning.psi = 2)

    summary(alr_Silt_spec_b_psi2)

alr\_Clay\_lab

    vario_alr_Clay_lab_b <- sample.variogram(data.spec[data.spec$layer == "B" & data.spec$set != "validation","alr_Clay_lab"], 
                                             locations = data.spec[data.spec$layer == "B"  & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                             lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Clay_lab_b[,c("lag.dist")], vario_alr_Clay_lab_b[,c("gamma")])

    alr_Clay_lab_b_psi2 <- georob(alr_Clay_lab ~ 1, data = data.spec[data.spec$layer == "B" & data.spec$set != "validation",], 
                                   locations = ~ POINT_X + POINT_Y, 
                                   variogram.model = "RMexp",
                                   param = c(variance = 1, nugget = 0.3, scale = 1500),
                                   fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                   tuning.psi = 2)

    summary(alr_Clay_lab_b_psi2)

alr\_Clay\_spec

    vario_alr_Clay_spec_b <- sample.variogram(data.spec[data.spec$layer == "B" & data.spec$set != "validation","alr_Clay_spec"], 
                                              locations = data.spec[data.spec$layer == "B" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                              lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_alr_Clay_spec_b[,c("lag.dist")], vario_alr_Clay_spec_b[,c("gamma")])

    alr_Clay_spec_b_psi2 <- georob(alr_Clay_spec ~ 1, data = data.spec[data.spec$layer == "B" & data.spec$set != "validation",], 
                                    locations = ~ POINT_X + POINT_Y, 
                                    variogram.model = "RMexp",
                                    param = c(variance = 1, nugget = 0.1, scale = 1500),
                                    fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                                    tuning.psi = 2)

    summary(alr_Clay_spec_b_psi2)

Ca lab

    vario_Ca_lab_b <- sample.variogram(data.spec[data.spec$layer == "B" & data.spec$set != "validation","Ca"], 
                                       locations = data.spec[data.spec$layer == "B" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                       lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_Ca_lab_b[,c("lag.dist")], vario_Ca_lab_b[,c("gamma")])

    Ca_lab_b_psi2 <- georob(Ca ~ 1, data = data.spec[data.spec$layer == "B" & data.spec$set != "validation",], 
                             locations = ~ POINT_X + POINT_Y, 
                             variogram.model = "RMexp",
                             param = c(variance = 200, nugget = 0.1, scale = 2000),
                             fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                             tuning.psi = 2)

    summary(Ca_lab_b_psi2)

specCa

    vario_specCa_b <- sample.variogram(data.spec[data.spec$layer == "B" & data.spec$set != "validation","specCa"], 
                                       locations = data.spec[data.spec$layer == "B" & data.spec$set != "validation", c("POINT_X", "POINT_Y")], 
                                       lag.dist.def = seq(0, 1500, by = 100), estimator = "matheron")
    plot(vario_specCa_b[,c("lag.dist")], vario_specCa_b[,c("gamma")])

    specCa_b_psi2 <- georob(specCa ~ 1, data = data.spec[data.spec$layer == "B" & data.spec$set != "validation",], 
                             locations = ~ POINT_X + POINT_Y, 
                             variogram.model = "RMexp",
                             param = c(variance = 150, nugget = 0.1, scale = 2000),
                             fit.param = default.fit.param(scale = TRUE, alpha = TRUE),
                             tuning.psi = 2)

    summary(specCa_b_psi2)

---- Validations layer A ----
=============================

Texture lab

    pred_alr_Silt_lab_a <- predict(alr_Silt_lab_a_psi2, 
                                   newdata = data.spec[data.spec$layer == "A" & data.spec$set == "validation",], 
                                   control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
    str(pred_alr_Silt_lab_a)


    pred_alr_Clay_lab_a <- predict(alr_Clay_lab_a_psi2, 
                                   newdata = data.spec[data.spec$layer == "A" & data.spec$set == "validation",], 
                                   control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
    str(pred_alr_Clay_lab_a)

Back-transformation to the fractions and precision measures
===========================================================

    dvr.Silt <- exp(pred_alr_Silt_lab_a$pred$pred + (0.5 * (pred_alr_Silt_lab_a$pred$var.target - pred_alr_Silt_lab_a$pred$cov.pred.target)))
    dvr.Clay <- exp(pred_alr_Clay_lab_a$pred$pred + (0.5 * (pred_alr_Clay_lab_a$pred$var.target - pred_alr_Clay_lab_a$pred$cov.pred.target)))
    dvn.pred <- 1 + dvr.Silt + dvr.Clay


    pred_Silt_lab_a <- 100 * (dvr.Silt/dvn.pred)
    pred_Clay_lab_a <- 100 * (dvr.Clay/dvn.pred)
    pred_Sand_lab_a <- 100/dvn.pred


    valR2_Silt_a <- (cor(pred_Silt_lab_a, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Silt"]))^2
    valrmse_Silt_a <- mean((pred_Silt_lab_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Silt"])^2)^0.5
    valme_Silt_a <- mean(pred_Silt_lab_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Silt"])



    valR2_Clay_a <- (cor(pred_Clay_lab_a, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Clay"]))^2
    valrmse_Clay_a <- mean((pred_Clay_lab_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Clay"])^2)^0.5
    valme_Clay_a <- mean(pred_Clay_lab_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Clay"])

    valR2_Sand_a <- (cor(pred_Sand_lab_a, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Sand"]))^2
    valrmse_Sand_a <- mean((pred_Sand_lab_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Sand"])^2)^0.5
    valme_Sand_a <- mean(pred_Sand_lab_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Sand"])

Texture spec

    pred_alr_Silt_spec_a <- predict(alr_Silt_spec_a_psi2, 
                                    newdata = data.spec[data.spec$layer == "A" & data.spec$set == "validation",], 
                                    control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE),
                                    param=c(variance = as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[2] + var.res.alr_silt_a)))

    str(pred_alr_Silt_spec_a)


    pred_alr_Clay_spec_a <- predict(alr_Clay_spec_a_psi2, 
                                    newdata = data.spec[data.spec$layer == "A" & data.spec$set == "validation",], 
                                    control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE),
                                    param=c(variance = as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[2] + var.res.alr_clay_a)))
    str(pred_alr_Clay_spec_a)

Back-transformation to the fractions and precision measures
===========================================================

    dvr.Silt_spec <- exp(pred_alr_Silt_spec_a$pred$pred + (0.5 * (pred_alr_Silt_spec_a$pred$var.target - pred_alr_Silt_spec_a$pred$cov.pred.target)))
    dvr.Clay_spec <- exp(pred_alr_Clay_spec_a$pred$pred + (0.5 * (pred_alr_Clay_spec_a$pred$var.target - pred_alr_Clay_spec_a$pred$cov.pred.target)))
    dvn.pred_spec <- 1 + dvr.Silt_spec + dvr.Clay_spec

    pred_Silt_spec_a <- 100 * (dvr.Silt_spec/dvn.pred_spec)
    pred_Clay_spec_a <- 100 * (dvr.Clay_spec/dvn.pred_spec)
    pred_Sand_spec_a <- 100/dvn.pred_spec


    valR2_spec_Silt_a <- (cor(pred_Silt_spec_a, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Silt"]))^2
    valrmse_spec_Silt_a <- mean((pred_Silt_spec_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Silt"])^2)^0.5
    valme_spec_Silt_a <- mean(pred_Silt_spec_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Silt"])

    valR2_spec_Clay_a <- (cor(pred_Clay_spec_a, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Clay"]))^2
    valrmse_spec_Clay_a <- mean((pred_Clay_spec_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Clay"])^2)^0.5
    valme_spec_Clay_a <- mean(pred_Clay_spec_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Clay"])

    valR2_spec_Sand_a <- (cor(pred_Sand_spec_a, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Sand"]))^2
    valrmse_spec_Sand_a <- mean((pred_Sand_spec_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Sand"])^2)^0.5
    valme_spec_Sand_a <- mean(pred_Sand_spec_a - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Sand"])

Ca lab

    pred_Ca_lab_a <- predict(Ca_lab_a_psi2, 
                             newdata = data.spec[data.spec$layer == "A" & data.spec$set == "validation",], 
                             control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
    str(pred_Ca_lab_a)

    valR2_lab_Ca_a <- (cor(pred_Ca_lab_a$pred$pred, data.spec[data.spec$layer == "A" & data.spec$set == "validation","Ca"]))^2
    valrmse_lab_Ca_a <- mean((pred_Ca_lab_a$pred$pred - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Ca"])^2)^0.5
    valme_lab_Ca_a <- mean(pred_Ca_lab_a$pred$pred - data.spec[data.spec$layer == "A" & data.spec$set == "validation","Ca"])

Ca spec

    pred_Ca_spec_a <- predict(specCa_a_psi2, 
        newdata = data.spec[data.spec$layer == 
            "A" & data.spec$set == "validation", 
            ], control = control.predict.georob(extended.output = TRUE, 
            full.covmat = TRUE), param = c(variance = as.numeric(specCa_a_psi2$variogram.object[[1]]$param[1]), 
            nugget = as.numeric(specCa_a_psi2$variogram.object[[1]]$param[3]), 
            scale = as.numeric(specCa_a_psi2$variogram.object[[1]]$param[4]), 
            snugget = as.numeric(specCa_a_psi2$variogram.object[[1]]$param[2] + 
                var.res.Ca_a)))
    str(pred_Ca_spec_a)

    valR2_spec_Ca_a <- (cor(pred_Ca_spec_a$pred$pred, 
        data.spec[data.spec$layer == "A" & data.spec$set == 
            "validation", "Ca"]))^2
    valrmse_spec_Ca_a <- mean((pred_Ca_spec_a$pred$pred - 
        data.spec[data.spec$layer == "A" & data.spec$set == 
            "validation", "Ca"])^2)^0.5
    valme_spec_Ca_a <- mean(pred_Ca_spec_a$pred$pred - 
        data.spec[data.spec$layer == "A" & data.spec$set == 
            "validation", "Ca"])

---- Validations layer B ----
=============================

Texture lab

    pred_alr_Silt_lab_b <- predict(alr_Silt_lab_b_psi2, 
                                   newdata = data.spec[data.spec$layer == "B" & data.spec$set == "validation",], 
                                   control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
    str(pred_alr_Silt_lab_b)

    pred_alr_Clay_lab_b <- predict(alr_Clay_lab_b_psi2, 
                                   newdata = data.spec[data.spec$layer == "B" & data.spec$set == "validation",], 
                                   control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
    str(pred_alr_Clay_lab_b)

Back-transformation to the fractions and precision measures
===========================================================

    dvr.Silt_b <- exp(pred_alr_Silt_lab_b$pred$pred + (0.5 * (pred_alr_Silt_lab_b$pred$var.target - pred_alr_Silt_lab_b$pred$cov.pred.target)))
    dvr.Clay_b <- exp(pred_alr_Clay_lab_b$pred$pred + (0.5 * (pred_alr_Clay_lab_b$pred$var.target - pred_alr_Clay_lab_b$pred$cov.pred.target)))
    dvn.pred_b <- 1 + dvr.Silt_b + dvr.Clay_b


    pred_Silt_lab_b <- 100 * (dvr.Silt_b/dvn.pred_b)
    pred_Clay_lab_b <- 100 * (dvr.Clay_b/dvn.pred_b)
    pred_Sand_lab_b <- 100/dvn.pred_b


    valR2_Silt_b <- (cor(pred_Silt_lab_b, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Silt"]))^2
    valrmse_Silt_b <- mean((pred_Silt_lab_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Silt"])^2)^0.5
    valme_Silt_b <- mean(pred_Silt_lab_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Silt"])

    valR2_Clay_b <- (cor(pred_Clay_lab_b, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Clay"]))^2
    valrmse_Clay_b <- mean((pred_Clay_lab_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Clay"])^2)^0.5
    valme_Clay_b <- mean(pred_Clay_lab_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Clay"])

    valR2_Sand_b <- (cor(pred_Sand_lab_b, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Sand"]))^2
    valrmse_Sand_b <- mean((pred_Sand_lab_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Sand"])^2)^0.5
    valme_Sand_b <- mean(pred_Sand_lab_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Sand"])

Texture spec

    pred_alr_Silt_spec_b <- predict(alr_Silt_spec_b_psi2, 
                                    newdata = data.spec[data.spec$layer == "B" & data.spec$set == "validation",], 
                                    control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE),
                              param=c(variance = as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[2] + var.res.alr_silt_b)))
    str(pred_alr_Silt_spec_b)



    pred_alr_Clay_spec_b <- predict(alr_Clay_spec_b_psi2, 
                                    newdata = data.spec[data.spec$layer == "B" & data.spec$set == "validation",], 
                                    control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE),
                              param=c(variance = as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[2] + var.res.alr_clay_b)))
    str(pred_alr_Clay_spec_b)

Back-transformation to the fractions and precision measures
===========================================================

    dvr.Silt_b_spec <- exp(pred_alr_Silt_spec_b$pred$pred + (0.5 * (pred_alr_Silt_spec_b$pred$var.target - pred_alr_Silt_spec_b$pred$cov.pred.target)))
    dvr.Clay_b_spec <- exp(pred_alr_Clay_spec_b$pred$pred + (0.5 * (pred_alr_Clay_spec_b$pred$var.target - pred_alr_Clay_spec_b$pred$cov.pred.target)))
    dvn.pred_b_spec <- 1 + dvr.Silt_b_spec + dvr.Clay_b_spec

    pred_Silt_spec_b <- 100 * (dvr.Silt_b_spec/dvn.pred_b_spec)
    pred_Clay_spec_b <- 100 * (dvr.Clay_b_spec/dvn.pred_b_spec)
    pred_Sand_spec_b <- 100/dvn.pred_b_spec

    valR2_spec_Silt_b <- (cor(pred_Silt_spec_b, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Silt"]))^2
    valrmse_spec_Silt_b <- mean((pred_Silt_spec_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Silt"])^2)^0.5
    valme_spec_Silt_b <- mean(pred_Silt_spec_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Silt"])

    valR2_spec_Clay_b <- (cor(pred_Clay_spec_b, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Clay"]))^2
    valrmse_spec_Clay_b <- mean((pred_Clay_spec_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Clay"])^2)^0.5
    valme_spec_Clay_b <- mean(pred_Clay_spec_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Clay"])

    valR2_spec_Sand_b <- (cor(pred_Sand_spec_b, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Sand"]))^2
    valrmse_spec_Sand_b <- mean((pred_Sand_spec_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Sand"])^2)^0.5
    valme_spec_Sand_b <- mean(pred_Sand_spec_b - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Sand"])

Ca lab

    pred_Ca_lab_b <- predict(Ca_lab_b_psi2, 
                             newdata = data.spec[data.spec$layer == "B" & data.spec$set == "validation",], 
                             control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE))
    str(pred_Ca_lab_b)

    valR2_lab_Ca_b <- (cor(pred_Ca_lab_b$pred$pred, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Ca"]))^2
    valrmse_lab_Ca_b <- mean((pred_Ca_lab_b$pred$pred - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Ca"])^2)^0.5
    valme_lab_Ca_b <- mean(pred_Ca_lab_b$pred$pred - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Ca"])

Ca spec

    pred_Ca_spec_b <- predict(specCa_b_psi2, 
                              newdata = data.spec[data.spec$layer == "B" & data.spec$set == "validation",], 
                              control = control.predict.georob(extended.output = TRUE, full.covmat = TRUE),
                              param=c(variance = as.numeric(specCa_b_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(specCa_b_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(specCa_b_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(specCa_b_psi2$variogram.object[[1]]$param[2] + var.res.Ca_b)))
    str(pred_Ca_spec_b)

    valR2_spec_Ca_b <- (cor(pred_Ca_spec_b$pred$pred, data.spec[data.spec$layer == "B" & data.spec$set == "validation","Ca"]))^2
    valrmse_spec_Ca_b <- mean((pred_Ca_spec_b$pred$pred - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Ca"])^2)^0.5
    valme_spec_Ca_b <- mean(pred_Ca_spec_b$pred$pred - data.spec[data.spec$layer == "B" & data.spec$set == "validation","Ca"])


    results_spatial_modeling_lab <- data.frame(Property = c("Sand", "Silt", "Clay", "Ca++", "Sand", "Silt", "Clay", "Ca++"),
                                               R2 = rep(NA, 8), 
                                               RMSE = rep(NA, 8),
                                               ME = rep(NA, 8),
                                               Layer = rep(c("A", "B"), each = 4))

    results_spatial_modeling_spec <- results_spatial_modeling_lab  

    results_spatial_modeling_lab$R2 <- c(valR2_Sand_a,
                                         valR2_Silt_a,
                                         valR2_Clay_a, 
                                         valR2_lab_Ca_a,
                                         valR2_Sand_b,
                                         valR2_Silt_b,
                                         valR2_Clay_b, 
                                         valR2_lab_Ca_b)
    results_spatial_modeling_lab$RMSE <- c(valrmse_Sand_a,
                                           valrmse_Silt_a,
                                           valrmse_Clay_a, 
                                           valrmse_lab_Ca_a,
                                           valrmse_Sand_b,
                                           valrmse_Silt_b,
                                           valrmse_Clay_b, 
                                           valrmse_lab_Ca_b)

    results_spatial_modeling_lab$ME <- c(valme_Sand_a,
                                           valme_Silt_a,
                                           valme_Clay_a, 
                                           valme_lab_Ca_a,
                                           valme_Sand_b,
                                           valme_Silt_b,
                                           valme_Clay_b, 
                                           valme_lab_Ca_b)

    results_spatial_modeling_spec$R2 <- c(valR2_spec_Sand_a,
                                          valR2_spec_Silt_a,
                                          valR2_spec_Clay_a, 
                                          valR2_spec_Ca_a,
                                          valR2_spec_Sand_b,
                                          valR2_spec_Silt_b,
                                          valR2_spec_Clay_b, 
                                          valR2_spec_Ca_b)

    results_spatial_modeling_spec$RMSE <- c(valrmse_spec_Sand_a,
                                            valrmse_spec_Silt_a,
                                            valrmse_spec_Clay_a, 
                                            valrmse_spec_Ca_a,
                                            valrmse_spec_Sand_b,
                                            valrmse_spec_Silt_b,
                                            valrmse_spec_Clay_b, 
                                            valrmse_spec_Ca_b)

    results_spatial_modeling_spec$ME <- c(valme_spec_Sand_a,
                                            valme_spec_Silt_a,
                                            valme_spec_Clay_a, 
                                            valme_spec_Ca_a,
                                            valme_spec_Sand_b,
                                            valme_spec_Silt_b,
                                            valme_spec_Clay_b, 
                                            valme_spec_Ca_b)

    results_parameters_weights_lab <- data.frame(Property = c("log(Silt/Sand)", "log(Clay/Sand","Ca++", "log(Silt/Sand)", "log(Clay/Sand", "Ca++"),
                                               Outliers = rep(NA, 6), 
                                               Layer = rep(c("A", "B"), each = 3))


    results_parameters_weights_spec <- results_parameters_weights_lab  

    results_parameters_weights_lab$Outliers <- c(sum(alr_Silt_lab_a_psi2$rweights<0.8),
                                               sum(alr_Clay_lab_a_psi2$rweights<0.8),
                                               sum(Ca_lab_a_psi2$rweights<0.8),
                                               sum(alr_Silt_lab_b_psi2$rweights<0.8),
                                               sum(alr_Clay_lab_b_psi2$rweights<0.8),
                                               sum(Ca_lab_b_psi2$rweights<0.8))

    results_parameters_weights_spec$Outliers <- c(sum(alr_Silt_spec_a_psi2$rweights<0.8),
                                                sum(alr_Clay_spec_a_psi2$rweights<0.8),
                                                sum(specCa_a_psi2$rweights<0.8),
                                                sum(alr_Silt_spec_b_psi2$rweights<0.8),
                                                sum(alr_Clay_spec_b_psi2$rweights<0.8),
                                                sum(specCa_b_psi2$rweights<0.8))

    results_parameters_weights <- cbind(results_parameters_weights_lab, results_parameters_weights_spec)
    results_spatial_modeling <- cbind(results_spatial_modeling_lab,results_spatial_modeling_spec) 

    write.table(results_parameters_weights, file = "results_parameters_weights.txt", row.names = FALSE, sep = "\t")
    write.table(results_spatial_modeling, file = "results_val_spatial_modeling.txt", row.names = FALSE, sep = "\t")

Plot the variogram for the layer A and B

    ### Get dataset concerning empirical semivariograms
    dataS  <- c( 'lab', 'spec' )
    layerS <- c( 'a', 'b' )
    propS  <- c( 'Silt', 'Clay', 'Ca' )


    # Prepare dataset
    for( i in 1:length(dataS)){  # dataset
      for( j in 1:length(layerS)){  # layer
        for( k in 1:length(propS)){  # property
         if( propS[k] != 'Ca'){
          if(i==1 & j==1 & k==1){
            dataBind  <-  as.data.frame( eval( parse(text=paste('vario_alr', propS[k],  dataS[i], layerS[j],  sep='_')) )[c('lag.dist', 'gamma')] )
            dataBind['PropS']   <-rep(propS[k],  dim(dataBind)[1]) 
            dataBind['dataS']   <-rep(dataS[i],  dim(dataBind)[1]) 
            dataBind['layerS']  <-rep(layerS[j],  dim(dataBind)[1])
                  
          }else{
            dataBindT <-  as.data.frame( eval( parse(text=paste('vario_alr', propS[k],  dataS[i], layerS[j],  sep='_')) )[c('lag.dist', 'gamma')] )
            dataBindT['PropS']   <-rep(propS[k],   dim(dataBindT)[1]) 
            dataBindT['dataS']   <-rep(dataS[i],   dim(dataBindT)[1]) 
            dataBindT['layerS']  <-rep(layerS[j],  dim(dataBindT)[1])

            dataBind  <- rbind(dataBind, dataBindT)
        
          }
         }else{
           if( dataS[i] == 'lab' ){
            dataBindT <- as.data.frame( eval( parse(text=paste('vario', propS[k],  dataS[i], layerS[j],  sep='_')) )[c('lag.dist', 'gamma')] )
            dataBindT['PropS']   <-rep(propS[k],   dim(dataBindT)[1]) 
            dataBindT['dataS']   <-rep(dataS[i],   dim(dataBindT)[1]) 
            dataBindT['layerS']  <-rep(layerS[j],  dim(dataBindT)[1])
            
            dataBind  <- rbind(dataBind, dataBindT)       

           }else{
            dataBindT  <- as.data.frame( eval( parse(text=paste('vario', '_', dataS[i], propS[k], '_', layerS[j],  sep='')) )[c('lag.dist', 'gamma')] )
            dataBindT['PropS']   <-rep(propS[k],   dim(dataBindT)[1]) 
            dataBindT['dataS']   <-rep(dataS[i],   dim(dataBindT)[1]) 
            dataBindT['layerS']  <-rep(layerS[j],  dim(dataBindT)[1])        
            
            dataBind  <- rbind(dataBind,  dataBindT)
         }
        }
       }
      }
     }
    colnames(dataBind) <- c("dist", "gamma", "PropS", "dataS", "layerS")

    ### Get dataset concerning fitted model
    for( i in 1:length(dataS)){  # dataset
      for( j in 1:length(layerS)){  # layer
        for( k in 1:length(propS)){  # property
          if( propS[k] != 'Ca'){
            if(i==1 & j==1 & k==1){
              georobT <- eval( parse(text=paste('alr', propS[k],  dataS[i], layerS[j], 'psi2',  sep='_')) )
              
              dataFitted <- variogramLine( vgm(psill = georobT$variogram.object[[1]]$param[1], 
                                                  model = "Exp", range = georobT$variogram.object[[1]]$param[4],
                                                  nugget = georobT$variogram.object[[1]]$param[3]), 2000 )
              dataFitted['PropS']   <-rep(propS[k],  dim(dataFitted)[1]) 
              dataFitted['dataS']   <-rep(dataS[i],  dim(dataFitted)[1]) 
              dataFitted['layerS']  <-rep(layerS[j], dim(dataFitted)[1])
              }else{
              georobT <- eval( parse(text=paste('alr', propS[k],  dataS[i], layerS[j], 'psi2',  sep='_')) )
              dataFittedT <- variogramLine( vgm(psill = georobT$variogram.object[[1]]$param[1], 
                                               model = "Exp", range = georobT$variogram.object[[1]]$param[4],
                                               nugget = georobT$variogram.object[[1]]$param[3]), 2000 )
              dataFittedT['PropS']   <-rep(propS[k],   dim(dataFittedT)[1]) 
              dataFittedT['dataS']   <-rep(dataS[i],   dim(dataFittedT)[1]) 
              dataFittedT['layerS']  <-rep(layerS[j],  dim(dataFittedT)[1])
              dataFitted  <- rbind(dataFitted, dataFittedT)
              
            }
          }else{
            if( dataS[i] == 'lab' ){
              georobT <- eval( parse(text=paste(propS[k], dataS[i], layerS[j], 'psi2',  sep='_')) )
              
              dataFittedT <- variogramLine( vgm(psill = georobT$variogram.object[[1]]$param[1], 
                                                model = "Exp", range = georobT$variogram.object[[1]]$param[4],
                                                nugget = georobT$variogram.object[[1]]$param[3]), 2000 )
              dataFittedT['PropS']   <-rep(propS[k],   dim(dataFittedT)[1]) 
              dataFittedT['dataS']   <-rep(dataS[i],   dim(dataBindT)[1]) 
              dataFittedT['layerS']  <-rep(layerS[j],  dim(dataFittedT)[1])
              dataFitted  <- rbind(dataFitted, dataFittedT)       
            }else{
              georobT  <- eval( parse(text=paste(dataS[i], propS[k], '_', layerS[j], '_', 'psi2',   sep='')) )
              
              dataFittedT <- variogramLine( vgm(psill = georobT$variogram.object[[1]]$param[1], 
                                                model = "Exp", range = georobT$variogram.object[[1]]$param[4],
                                                nugget = georobT$variogram.object[[1]]$param[3]), 2000 )
              
              dataFittedT['PropS']   <-rep(propS[k],   dim(dataFittedT)[1]) 
              dataFittedT['dataS']   <-rep(dataS[i],   dim(dataFittedT)[1]) 
              dataFittedT['layerS']  <-rep(layerS[j],  dim(dataFittedT)[1])        
              
              dataFitted  <- rbind(dataFitted,  dataFittedT)
              
            }
          }
        }
      }
    }

    ### merge empirical and fitted values before plotting
    dataBindAll  <-  dataBind
    dataBindAll['target'] <-  rep( 'empir', dim(dataBindAll)[1] )

    dataBindAllT  <-  dataFitted
    dataBindAllT['target'] <-  rep( 'fitt', dim(dataFitted)[1] )

    dataBindAll <- rbind(dataBindAll, dataBindAllT)

    ### Model Parameters
    dataModel <- read.table( "ModelParam.txt",  sep = "\t", header=T)

    ### Plot without facets
    dataS2 <- unique( as.vector(dataBindAll[['dataS']]) )
    PropS2 <- unique( as.vector(dataBindAll[['PropS']]) )
    PropS2P <- rev( c("'Ca'^'++'", 'log(Clay / Sand)', 'log(Silt / Sand)' ) )

    # get limits for y-axis in each plot
    nbreaks = 4
    limAll = list()
    for( j in 1:length(layerS)){  # layer
        
     limAllT = list()
     for( k in 1:length(PropS2)){  # property
      for( i in 1:length(dataS2)){  # dataset
          dataSelec <- scales::pretty_breaks(n = nbreaks)( dataBindAll[ ( dataBindAll[['dataS']] == dataS2[i] & dataBindAll[['layerS']] == layerS[j] &
                                                                          dataBindAll[['PropS']] == PropS2[k] ), ][['gamma']] )
          limAllT[[ (k*2)-(2-i) ]] <- max(dataSelec)
        }
      }    
        
      limAll[j] = list(limAllT)
    }

    # get vector with coordinates for text
    limAllFF=list()
    for(j in 1:length(layerS)){ # layer
      limAllF = c()
      for( k in 1:length(PropS2)){  # property
        for( i in 1:length(dataS2)){  # dataset
          limAllF <-  c(limAllF, limAll[[j]][[(k*2)-(2-i)]])
        }
      }
      limAllFF[[j]] <- limAllF
    }

    # color and location
    col1 <- 'blue4'
    col2 <- 'orangered3'
    xdist <- 1400
    xdistM1 <- 1375
    xdistM2 <- Inf
    fontSize1 <- 8
    fontSize2 <- 8
    fontSize3 <- 2.8
    lineSize1 <- 0.6
    borderSize1 <- 0.3
    pointSize1 <- 1.5

    # Adjust the ploting option
    calcF <- function(data, paramP){ 
              if(paramP==1){
                val <- 0.375*data
               }else{
                if(paramP==2){
                  val <- 0.27*data
                 }else{
                  if(paramP==3){
                    val <- 0.15*data
                  }else{
                    if(paramP==4){
                     val <- 0.06*data
                    }else{
                      if(paramP==5)  
                       val <- 0.47*data
                      else{
                        val <- 0.50*data
                      }
                    }    
                  }
                }
              }
             return(val) 
        }
      
    #
    for(i in 1:length(layerS)){  # layer

          
          dataP      <- dataBindAll[ dataBindAll['layerS'] == layerS[i] & dataBindAll['dataS'] == dataS2[1] & dataBindAll['PropS'] == PropS2[1], ]
          datamodelP <- dataModel[ dataModel['layerS'] == layerS[i] & dataModel['dataS'] == dataS2[1] & dataModel['PropS'] == PropS2[1], ]
          
          #
          P1 <- 
            ggplot( dataP, aes(x=dist, y=gamma) )+
            geom_point( data=subset(dataP, target == "empir"), aes(colour=target),  size=pointSize1, shape=1)+   # , stro1e = 1
            geom_line(  data=subset(dataP, target == "fitt"),  aes(x=dist, y=gamma, colour=target), size=lineSize1, linetype="longdash")+
            scale_color_manual(values=c(col1, col2), labels = c("Empirical", "Fitted"))+
            labs(x='', y=paste('Semivariance of ', PropS2P[1]) )+
            annotate( 'text', x=xdist, y=calcF(limAllFF[[i]][(1*2)-(1-2)],5),  label="italic('Fitted variogram:')", hjust = 0, vjust = 1, parse=T, size=fontSize3)+
            annotate("rect", xmin =xdistM1, xmax =xdistM2, ymin =-Inf, ymax=calcF(limAllFF[[i]][(1*2)-(1-2)],6), colour='black', fill = NA, size=0.3)+
            geom_text(data=datamodelP, aes( label=paste('tau^2', '==', tau),     x=xdist,  y=calcF(limAllFF[[i]][(1*2)-(1-2)],1)), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2', '==', sigma), x=xdist,  y=calcF(limAllFF[[i]][(1*2)-(1-2)],2)), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('alpha', '==', alpha),   x=xdist,  y=calcF(limAllFF[[i]][(1*2)-(1-2)],3)), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2 / ( sigma^2 + tau^2 )', '==', sigmaTau), x=xdist, y=calcF(limAllFF[[i]][(1*2)-(1-2)],4)), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            ggtitle( paste('Laboratory-based | Layer ', c('A','B')[i], ' ', c('(0-0.2', '(0.8-1.0')[i], ', n = 345)', sep='') )+
            scale_x_continuous(breaks = pretty(dataP$dist, n = 4), limits=c(0,NA))+
            scale_y_continuous(breaks = pretty(dataP$gamma, n = 4), limits=c(0,NA))+
            theme( axis.title.x = element_text(size=fontSize2, colour="black", family="Arial"),
                   axis.title.y = element_text(size=fontSize2, colour="black", family="Arial"), 
                   axis.text.x  = element_text(size=fontSize1, colour="black", family="Arial"), 
                   axis.text.y  = element_text(size=fontSize1, colour="black", family="Arial"),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black",  fill="transparent", size=borderSize1, linetype='solid'),
                   plot.background =  element_blank(),
                   plot.title = element_text(colour = "black", size=fontSize2, hjust = 0.5, vjust = 0, family="Arial", face='plain'),
                   axis.line    = element_blank(), 
                   axis.ticks.y=element_line(size = 0.3),
                   axis.ticks.x = element_line(size = 0.3),
                   legend.position ="none",
                   plot.margin= unit(c(0,0.1,0,0), "cm"))

          dataP      <- dataBindAll[ dataBindAll['layerS'] == layerS[i] & dataBindAll['dataS'] == dataS2[2] & dataBindAll['PropS'] == PropS2[1], ]
          datamodelP <- dataModel[ dataModel['layerS'] == layerS[i] & dataModel['dataS'] == dataS2[2] & dataModel['PropS'] == PropS2[1], ]
          

          P2 <- 
            ggplot( dataP, aes(x=dist, y=gamma) )+
            geom_point( data=subset(dataP, target == "empir"), aes(colour=target),  size=pointSize1, shape=1)+   # , stroke = 1
            geom_line(  data=subset(dataP, target == "fitt"),  aes(x=dist, y=gamma, colour=target), size=lineSize1, linetype="dashed")+
            scale_color_manual(values=c(col1, col2), labels = c("Empirical", "Fitted"))+
            labs(x='', y=paste('Semivariance of ', PropS2P[1]))+
            annotate( 'text', x=xdist, y=calcF(limAllFF[[i]][(1*2)-(2-2)],5),  label="italic('Fitted variogram:')", hjust = 0, vjust = 1, parse=T, size=fontSize3)+
            annotate("rect", xmin =xdistM1, xmax =xdistM2, ymin =-Inf, ymax=calcF(limAllFF[[i]][(1*2)-(2-2)],6), colour='black', fill = NA, size=0.3)+
            geom_text(data=datamodelP, aes( label=paste('tau^2', '==', tau),     x=xdist,  y=calcF(limAllFF[[i]][(1*2)-(2-2)],1) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2', '==', sigma), x=xdist,  y=calcF(limAllFF[[i]][(1*2)-(2-2)],2) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('alpha', '==', alpha),   x=xdist,  y=calcF(limAllFF[[i]][(1*2)-(2-2)],3) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2 / ( sigma^2 + tau^2 )', '==', sigmaTau), x=xdist, y=calcF(limAllFF[[i]][(1*2)-(2-2)],4) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            ggtitle( paste('Vis-NIR-aided | Layer ', c('A','B')[i], ' ', c('(0-0.2', '(0.8-1.0')[i], ', n = 345)', sep='') )+
            scale_x_continuous(breaks = pretty(dataP$dist, n = 4), limits=c(0,NA))+
            scale_y_continuous(breaks = pretty(dataP$gamma, n = 4), limits=c(0,NA))+
            theme( axis.title.x = element_text(size=fontSize2, colour="black", family="Arial"),
                   axis.title.y = element_text(size=fontSize2, colour="black", family="Arial"), 
                   axis.text.x  = element_text(size=fontSize1, colour="black", family="Arial"), 
                   axis.text.y  = element_text(size=fontSize1, colour="black", family="Arial"),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black",  fill="transparent", size=borderSize1, linetype='solid'),
                   plot.background =  element_blank(),
                   plot.title = element_text(colour = "black", size=fontSize2, hjust = 0.5, vjust = 0, family="Arial", face='plain'),
                   axis.line    = element_blank(),
                   axis.ticks.y=element_line(size = 0.3),
                   axis.ticks.x = element_line(size = 0.3),
                   legend.position ="none",
                   plot.margin= unit(c(0,0.1,0,0), "cm"))
          #      
          dataP      <- dataBindAll[ dataBindAll['layerS'] == layerS[i] & dataBindAll['dataS'] == dataS2[1] & dataBindAll['PropS'] == PropS2[2], ]
          datamodelP <- dataModel[ dataModel['layerS'] == layerS[i] & dataModel['dataS'] == dataS2[1] & dataModel['PropS'] == PropS2[2], ]
          
          P3 <- 
            ggplot( dataP, aes(x=dist, y=gamma) )+
            geom_point( data=subset(dataP, target == "empir"), aes(colour=target),  size=pointSize1, shape=1)+   # , stroke = 1
            geom_line(  data=subset(dataP, target == "fitt"),  aes(x=dist, y=gamma, colour=target), size=lineSize1, linetype="dashed")+
            scale_color_manual(values=c(col1, col2), labels = c("Empirical", "Fitted"))+
            labs(x='', y=paste('Semivariance of ', PropS2P[2]))+
            annotate( 'text', x=xdist, y=calcF(limAllFF[[i]][(2*2)-(2-1)],5),  label="italic('Fitted variogram:')", size=fontSize3, hjust = 0, vjust = 1, parse=T)+
            annotate("rect", xmin =xdistM1, xmax =xdistM2, ymin =-Inf, ymax=calcF(limAllFF[[i]][(2*2)-(2-1)],6), colour='black', fill = NA, size=0.3)+
            geom_text(data=datamodelP, aes( label=paste('tau^2', '==', tau),     x=xdist,  y=calcF(limAllFF[[i]][(2*2)-(2-1)],1) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2', '==', sigma), x=xdist,  y=calcF(limAllFF[[i]][(2*2)-(2-1)],2) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('alpha', '==', alpha),   x=xdist,  y=calcF(limAllFF[[i]][(2*2)-(2-1)],3) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2 / ( sigma^2 + tau^2 )', '==', sigmaTau), x=xdist, y=calcF(limAllFF[[i]][(2*2)-(2-1)],4) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            ggtitle( paste('Laboratory-based | Layer ', c('A','B')[i], ' ', c('(0-0.2', '(0.8-1.0')[i], ', n = 345)', sep='') )+
            scale_x_continuous(breaks = pretty(dataP$dist,  n = 4), limits=c(0,NA))+
            scale_y_continuous(breaks = pretty(dataP$gamma, n = 4), limits=c(0,NA))+
            theme( axis.title.x = element_text(size=fontSize2, colour="black", family="Arial"),
                   axis.title.y = element_text(size=fontSize2, colour="black", family="Arial"), 
                   axis.text.x  = element_text(size=fontSize1, colour="black", family="Arial"), 
                   axis.text.y  = element_text(size=fontSize1, colour="black", family="Arial"),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black",  fill="transparent", size=borderSize1, linetype='solid'),
                   plot.background =  element_blank(),
                   plot.title = element_text(colour = "black", size=fontSize2, hjust = 0.5, vjust = 0, family="Arial", face='plain'),
                   axis.line    = element_blank(), 
                   axis.ticks.y=element_line(size = 0.3),
                   axis.ticks.x = element_line(size = 0.3),
                   legend.position ="none",
                   plot.margin= unit(c(0,0.1,0,0), "cm"))  

          #      
          dataP      <- dataBindAll[ dataBindAll['layerS'] == layerS[i] & dataBindAll['dataS'] == dataS2[2] & dataBindAll['PropS'] == PropS2[2], ]
          datamodelP <- dataModel[ dataModel['layerS'] == layerS[i] & dataModel['dataS'] == dataS2[2] & dataModel['PropS'] == PropS2[2], ]

          P4 <- 
            ggplot( dataP, aes(x=dist, y=gamma) )+
            geom_point( data=subset(dataP, target == "empir"), aes(colour=target),  size=pointSize1, shape=1)+   # , stroke = 1
            geom_line(  data=subset(dataP, target == "fitt"),  aes(x=dist, y=gamma, colour=target), size=lineSize1, linetype="dashed")+
            scale_color_manual(values=c(col1, col2), labels = c("Empirical", "Fitted"))+
            labs(x='', y=paste('Semivariance of ', PropS2P[2]))+
            annotate( 'text', x=xdist, y=calcF(limAllFF[[i]][(2*2)-(2-2)],5),  label="italic('Fitted variogram:')", size=fontSize3, hjust = 0, vjust = 1, parse=T)+
            annotate("rect", xmin =xdistM1, xmax =xdistM2, ymin =-Inf, ymax=calcF(limAllFF[[i]][(2*2)-(2-2)],6), colour='black', fill = NA, size=0.3)+
            geom_text(data=datamodelP, aes( label=paste('tau^2', '==', tau),     x=xdist,  y=calcF(limAllFF[[i]][(2*2)-(2-2)],1) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2', '==', sigma), x=xdist,  y=calcF(limAllFF[[i]][(2*2)-(2-2)],2) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('alpha', '==', alpha),   x=xdist,  y=calcF(limAllFF[[i]][(2*2)-(2-2)],3) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2 / ( sigma^2 + tau^2 )', '==', sigmaTau), x=xdist, y=calcF(limAllFF[[i]][(2*2)-(2-2)],4) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            ggtitle( paste('Vis-NIR-aided | Layer ', c('A','B')[i], ' ', c('(0-0.2', '(0.8-1.0')[i], ', n = 345)', sep='') )+
            scale_x_continuous(breaks = pretty(dataP$dist, n = 4), limits=c(0,NA))+
            scale_y_continuous(breaks = pretty(dataP$gamma, n = 4), limits=c(0,NA))+
            theme( title = element_text(size=10, colour="black"),
                   axis.title.x = element_text(size=fontSize2, colour="black", family="Arial"),
                   axis.title.y = element_text(size=fontSize2, colour="black", family="Arial"), 
                   axis.text.x  = element_text(size=fontSize1, colour="black", family="Arial"), 
                   axis.text.y  = element_text(size=fontSize1, colour="black", family="Arial"),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black",  fill="transparent", size=borderSize1, linetype='solid'),
                   plot.background =  element_blank(),
                   plot.title = element_text(colour = "black", size=fontSize2, hjust = 0.5, vjust = 0, family="Arial", face='plain'),
                   axis.line    = element_blank(), 
                   axis.ticks.y=element_line(size = 0.3),
                   axis.ticks.x = element_line(size = 0.3),
                   legend.position ="none",
                   plot.margin= unit(c(0,0.1,0,0), "cm")) 

          #      
          dataP      <- dataBindAll[ dataBindAll['layerS'] == layerS[i] & dataBindAll['dataS'] == dataS2[1] & dataBindAll['PropS'] == PropS2[3], ]
          datamodelP <- dataModel[ dataModel['layerS'] == layerS[i] & dataModel['dataS'] == dataS2[1] & dataModel['PropS'] == PropS2[3], ]

          P5 <- 
            ggplot( dataP, aes(x=dist, y=gamma) )+
            geom_point( data=subset(dataP, target == "empir"), aes(colour=target),  size=pointSize1, shape=1)+   # , stroke = 1
            geom_line(  data=subset(dataP, target == "fitt"),  aes(x=dist, y=gamma, colour=target), size=lineSize1, linetype="dashed")+
            scale_color_manual(values=c(col1, col2), labels = c("Empirical", "Fitted"))+
            labs(x="Distance [m]", y=expression(Semivariance~of~Ca^'++') )+
            annotate( 'text', x=xdist, y=calcF(limAllFF[[i]][(3*2)-(2-1)],5),  label="italic('Fitted variogram:')", size=fontSize3, hjust = 0, vjust = 1, parse=T)+
            annotate("rect", xmin =xdistM1, xmax =xdistM2, ymin =-Inf, ymax=calcF(limAllFF[[i]][(3*2)-(2-1)],6), colour='black', fill = NA, size=0.3)+
            geom_text(data=datamodelP, aes( label=paste('tau^2', '==', tau),     x=xdist,  y=calcF(limAllFF[[i]][(3*2)-(2-1)],1) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2', '==', sigma), x=xdist,  y=calcF(limAllFF[[i]][(3*2)-(2-1)],2) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('alpha', '==', alpha),   x=xdist,  y=calcF(limAllFF[[i]][(3*2)-(2-1)],3) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2 / ( sigma^2 + tau^2 )', '==', sigmaTau), x=xdist, y=calcF(limAllFF[[i]][(3*2)-(2-1)],4) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            ggtitle( paste('Laboratory-based | Layer ', c('A','B')[i], ' ', c('(0-0.2', '(0.8-1.0')[i], ', n = 345)', sep='') )+
            scale_x_continuous(breaks = pretty(dataP$dist,  n = 4), limits=c(0,NA))+
            scale_y_continuous(breaks = pretty(dataP$gamma, n = 4), limits=c(0,NA))+
            theme( axis.title.x = element_text(size=fontSize2, colour="black", family="Arial"),
                   axis.title.y = element_text(size=fontSize2, colour="black", family="Arial"), 
                   axis.text.x  = element_text(size=fontSize1, colour="black", family="Arial"), 
                   axis.text.y  = element_text(size=fontSize1, colour="black", family="Arial"),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black",  fill="transparent", size=borderSize1, linetype='solid'),
                   plot.background =  element_blank(),
                   plot.title = element_text(colour = "black", size=fontSize2, hjust = 0.5, vjust = 0, family="Arial", face='plain'),
                   axis.line    = element_blank(), 
                   axis.ticks.y=element_line(size = 0.3),
                   axis.ticks.x = element_line(size = 0.3),
                   legend.position ="none",
                   plot.margin= unit(c(0,0.1,0,0), "cm")) 
     
          #      
          dataP      <- dataBindAll[ dataBindAll['layerS'] == layerS[i] & dataBindAll['dataS'] == dataS2[2] & dataBindAll['PropS'] == PropS2[3], ]
          datamodelP <- dataModel[ dataModel['layerS'] == layerS[i] & dataModel['dataS'] == dataS2[2] & dataModel['PropS'] == PropS2[3], ]

          P6 <- 
            ggplot( dataP, aes(x=dist, y=gamma) )+
            geom_point( data=subset(dataP, target == "empir"), aes(colour=target),  size=pointSize1, shape=1)+   # , stroke = 1
            geom_line(  data=subset(dataP, target == "fitt"),  aes(x=dist, y=gamma, colour=target), size=lineSize1, linetype="dashed")+
            scale_color_manual(values=c(col1, col2), labels = c("Empirical", "Fitted"))+
            labs(x="Distance [m]", y=expression(Semivariance~of~Ca^'++') )+
            annotate( 'text', x=xdist, y=calcF(limAllFF[[i]][(3*2)-(2-2)],5),  label="italic('Fitted variogram:')", size=fontSize3, hjust = 0, vjust = 1, parse=T)+
            annotate("rect", xmin =xdistM1, xmax =xdistM2, ymin =-Inf, ymax=calcF(limAllFF[[i]][(3*2)-(2-2)],6), colour='black', fill = NA, size=0.3)+
            geom_text(data=datamodelP, aes( label=paste('tau^2', '==', tau),     x=xdist,  y=calcF(limAllFF[[i]][(3*2)-(2-2)],1) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2', '==', sigma), x=xdist,  y=calcF(limAllFF[[i]][(3*2)-(2-2)],2) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('alpha', '==', alpha),   x=xdist,  y=calcF(limAllFF[[i]][(3*2)-(2-2)],3) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            geom_text(data=datamodelP, aes( label=paste('sigma^2 / ( sigma^2 + tau^2 )', '==', sigmaTau), x=xdist, y=calcF(limAllFF[[i]][(3*2)-(2-2)],4) ), parse=T, hjust = 0, vjust = 1, size=fontSize3 )+
            ggtitle( paste('Vis-NIR-aided | Layer ', c('A','B')[i], ' ', c('(0-0.2', '(0.8-1.0')[i], ', n = 345)', sep='') )+
            scale_x_continuous(breaks = pretty(dataP$dist,  n = 4), limits=c(0,NA))+
            scale_y_continuous(breaks = pretty(dataP$gamma, n = 4), limits=c(0,NA))+
            theme( axis.title.x = element_text(size=fontSize2, colour="black", family="Arial"),
                   axis.title.y = element_text(size=fontSize2, colour="black", family="Arial"), 
                   axis.text.x  = element_text(size=fontSize1, colour="black", family="Arial"), 
                   axis.text.y  = element_text(size=fontSize1, colour="black", family="Arial"),
                   panel.background = element_blank(),
                   panel.border = element_rect(colour = "black",  fill="transparent", size=borderSize1, linetype='solid'),
                   plot.background =  element_blank(),
                   plot.title = element_text(colour = "black", size=fontSize2, hjust = 0.5, vjust = 0, family="Arial", face='plain'),
                   axis.line    = element_blank(), 
                   axis.ticks.y=element_line(size = 0.3),
                   axis.ticks.x = element_line(size = 0.3),
                   legend.position ="none",
                   plot.margin= unit(c(0,0.1,0,0), "cm"))
          
          
          pAll <-plot_grid( P1,P2,P3,P4,P5,P6, labels = LETTERS[1:6], ncol = 2, scale=rep(0.98,6), label_fontfamily='Arial', label_size=10)
          save_plot(paste('SemiV_L', i, '.pdf'),  pAll, ncol = 2, base_height = 6.5, base_width = 3.6, dpi=600)
    }

--- Final spatial predictions ----
==================================

    pol2raster <-function(r,  nrows = 10, ncols = 10){
      ext<-raster(extent(r), nrows, ncols)
      crs(ext)<-crs(r)
      fr <-rasterize(r,ext,field=1,update=T)
      return(fr)
    }

    shape <- readOGR(dsn = workingd, 
                     layer = "lim_poly")
    summary(shape)
    spplot(shape)

    rncol <- (extent(shape)@ymax - extent(shape)@ymin)/10
    rnrow <- (extent(shape)@xmax - extent(shape)@xmin)/10

    rasg <- pol2raster(shape, rncol, rnrow)


    rasg <- raster::resample(rasg, 
                             raster(resolution = 10, ext = extent(shape)), 
                             method="ngb")
    rasg
    plot(rasg)

    sppx <- as(rasg, "SpatialPixelsDataFrame")
    colnames(sppx@coords) <- c("POINT_X", "POINT_Y")

Layer A

    alr_Silt_lab_a_map <- predict(alr_Silt_lab_a_psi2, sppx)
    alr_Clay_lab_a_map <- predict(alr_Clay_lab_a_psi2, sppx)

    Silt_lab_a_map <- Clay_lab_a_map <- Sand_lab_a_map <- alr_Silt_lab_a_map

    Silt_lab_a_map$pred <- 100 * exp(alr_Silt_lab_a_map$pred)/(1 + exp(alr_Silt_lab_a_map$pred) + exp(alr_Clay_lab_a_map$pred))
    Clay_lab_a_map$pred <- 100 * exp(alr_Clay_lab_a_map$pred)/(1 + exp(alr_Silt_lab_a_map$pred) + exp(alr_Clay_lab_a_map$pred))
    Sand_lab_a_map$pred <- 100/(1 + exp(alr_Silt_lab_a_map$pred) + exp(alr_Clay_lab_a_map$pred))


    alr_Clay_spec_a_map <- predict(alr_Clay_spec_a_psi2, sppx,
                              param=c(variance = as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Clay_spec_a_psi2$variogram.object[[1]]$param[2] + var.res.alr_clay_a)))
    alr_Silt_spec_a_map <- predict(alr_Silt_spec_a_psi2, sppx,
                              param=c(variance = as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Silt_spec_a_psi2$variogram.object[[1]]$param[2] + var.res.alr_silt_a)))

    Silt_spec_a_map <- Clay_spec_a_map <- Sand_spec_a_map <- alr_Silt_lab_a_map

    Silt_spec_a_map$pred <- 100 * exp(alr_Silt_spec_a_map$pred)/(1 + exp(alr_Silt_spec_a_map$pred) + exp(alr_Clay_spec_a_map$pred))
    Clay_spec_a_map$pred <- 100 * exp(alr_Clay_spec_a_map$pred)/(1 + exp(alr_Silt_spec_a_map$pred) + exp(alr_Clay_spec_a_map$pred))
    Sand_spec_a_map$pred <- 100/(1 + exp(alr_Silt_spec_a_map$pred) + exp(alr_Clay_spec_a_map$pred))


    Ca_lab_a_map <- predict(Ca_lab_a_psi2, sppx)
    Ca_spec_a_map <- predict(specCa_a_psi2, sppx,
                              param=c(variance = as.numeric(specCa_a_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(specCa_a_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(specCa_a_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(specCa_a_psi2$variogram.object[[1]]$param[2] + var.res.Ca_a)))

Layer B

    alr_Silt_lab_b_map <- predict(alr_Silt_lab_b_psi2, sppx)
    alr_Silt_spec_b_map <- predict(alr_Silt_spec_b_psi2, sppx,
                              param=c(variance = as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Silt_spec_b_psi2$variogram.object[[1]]$param[2] + var.res.alr_silt_b)))

    alr_Clay_lab_b_map <- predict(alr_Clay_lab_b_psi2, sppx)
    alr_Clay_spec_b_map <- predict(alr_Clay_spec_b_psi2, sppx,
                              param=c(variance = as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(alr_Clay_spec_b_psi2$variogram.object[[1]]$param[2] + var.res.alr_clay_b)))


    Silt_lab_b_map <- Clay_lab_b_map <- Sand_lab_b_map <- alr_Silt_lab_b_map

    Silt_lab_b_map$pred <- 100 * exp(alr_Silt_lab_b_map$pred)/(1 + exp(alr_Silt_lab_b_map$pred) + exp(alr_Clay_lab_b_map$pred))
    Clay_lab_b_map$pred <- 100 * exp(alr_Clay_lab_b_map$pred)/(1 + exp(alr_Silt_lab_b_map$pred) + exp(alr_Clay_lab_b_map$pred))
    Sand_lab_b_map$pred <- 100/(1 + exp(alr_Silt_lab_b_map$pred) + exp(alr_Clay_lab_b_map$pred))


    Silt_spec_b_map <- Clay_spec_b_map <- Sand_spec_b_map <- alr_Silt_lab_b_map

    Silt_spec_b_map$pred <- 100 * exp(alr_Silt_spec_b_map$pred)/(1 + exp(alr_Silt_spec_b_map$pred) + exp(alr_Clay_spec_b_map$pred))
    Clay_spec_b_map$pred <- 100 * exp(alr_Clay_spec_b_map$pred)/(1 + exp(alr_Silt_spec_b_map$pred) + exp(alr_Clay_spec_b_map$pred))
    Sand_spec_b_map$pred <- 100/(1 + exp(alr_Silt_spec_b_map$pred) + exp(alr_Clay_spec_b_map$pred))


    Ca_lab_b_map <- predict(Ca_lab_b_psi2, sppx)
    Ca_spec_b_map <- predict(specCa_b_psi2, sppx,
                              param=c(variance = as.numeric(specCa_b_psi2$variogram.object[[1]]$param[1]), 
                                      nugget =  as.numeric(specCa_b_psi2$variogram.object[[1]]$param[3]), 
                                      scale =  as.numeric(specCa_b_psi2$variogram.object[[1]]$param[4]),
                                      snugget =  as.numeric(specCa_b_psi2$variogram.object[[1]]$param[2] + var.res.Ca_b)))

Final Results

    pred_layer_A <- data.frame(rbind(coordinates(Silt_spec_a_map), coordinates(Silt_spec_a_map), coordinates(Silt_spec_a_map)),
                               layer = "Layer A (0 - 0.2 m)",
                               method = c(rep("Laboratory-based",nrow(Silt_spec_a_map)), 
                                          rep("vis-NIR augmented", nrow(Silt_spec_a_map)),
                                          rep("Differences between maps", nrow(Silt_spec_a_map))), 
                               Silt = c(Silt_lab_a_map$pred, Silt_spec_a_map$pred, Silt_lab_a_map$pred - Silt_spec_a_map$pred),
                               Clay = c(Clay_lab_a_map$pred, Clay_spec_a_map$pred, Clay_lab_a_map$pred - Clay_spec_a_map$pred),
                               Sand = c(Sand_lab_a_map$pred, Sand_spec_a_map$pred, Sand_lab_a_map$pred - Sand_spec_a_map$pred),
                               Ca = c(Ca_lab_a_map$pred, Ca_spec_a_map$pred, Ca_lab_a_map$pred - Ca_spec_a_map$pred),
                               row.names = 1:(3*nrow(Silt_spec_a_map)))



    pred_layer_B <- data.frame(rbind(coordinates(Silt_spec_b_map), coordinates(Silt_spec_b_map), coordinates(Silt_spec_b_map)),
                               layer = "Layer B (0.8 - 1.0 m)",
                               method = c(rep("Laboratory-based",nrow(Silt_spec_b_map)), 
                                          rep("vis-NIR augmented", nrow(Silt_spec_b_map)),
                                          rep("Differences between maps", nrow(Silt_spec_b_map))), 
                               Silt = c(Silt_lab_b_map$pred, Silt_spec_b_map$pred, Silt_lab_b_map$pred - Silt_spec_b_map$pred),
                               Clay = c(Clay_lab_b_map$pred, Clay_spec_b_map$pred, Clay_lab_b_map$pred - Clay_spec_b_map$pred),
                               Sand = c(Sand_lab_b_map$pred, Sand_spec_b_map$pred, Sand_lab_b_map$pred - Sand_spec_b_map$pred),
                               Ca = c(Ca_lab_b_map$pred, Ca_spec_b_map$pred, Ca_lab_b_map$pred - Ca_spec_b_map$pred),
                               row.names = 1:(3*nrow(Silt_spec_b_map)))

    preditions_final <- data.frame(rbind(pred_layer_A, pred_layer_B))

    pixpreds <- rbind(data.frame("Laboratorybased" = c(Sand_lab_a_map$pred, 
                                                        Silt_lab_a_map$pred, 
                                                        Clay_lab_a_map$pred,
                                                        Ca_lab_a_map$pred),
                                 "visNIRaided" = c(Sand_spec_a_map$pred, 
                                                     Silt_spec_a_map$pred, 
                                                     Clay_spec_a_map$pred,
                                                     Ca_spec_a_map$pred),
                                 Layer = "Layer A",
                                 Property = rep(c("Sand", "Silt", "Clay", "Ca"), each = length(Clay_spec_a_map$pred))),
                      data.frame("Laboratorybased" = c(Sand_lab_b_map$pred, 
                                                        Silt_lab_b_map$pred, 
                                                        Clay_lab_b_map$pred,
                                                        Ca_lab_b_map$pred),
                                 "visNIRaided" = c(Sand_spec_b_map$pred, 
                                                     Silt_spec_b_map$pred, 
                                                     Clay_spec_b_map$pred,
                                                     Ca_spec_b_map$pred),
                                 Layer = "Layer B",
                                 Property = rep(c("Sand", "Silt", "Clay", "Ca"), each = length(Clay_spec_b_map$pred))))

Plot results

    d2sand <- ggplot(pixpreds[pixpreds$Property == "Sand",], aes(Laboratorybased, visNIRaided)) + theme_bw() + 
      labs(y = "Sand content, % (vis-NIR augmented)", x = "Sand content, % (Laboratory-based)") +
      geom_point (alpha = 0.5, size = 3) + geom_abline(intercept = 0, slope = 1, colour = "red") +
      xlim(10, 90) + ylim(10, 90) +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.y  = element_text(angle=90, vjust =0.5, hjust =0.5, size=14), legend.position="top", legend.title=NULL) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.x  = element_text(angle=0, vjust=0, size=14)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=18),
            legend.title = element_text(face="bold", colour = grey(0.2), size=18)) +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold",size = 16, colour = "black", angle = 0)) 

    d2sand + facet_wrap(~ Layer, nrow = 1)


    d2silt <- ggplot(pixpreds[pixpreds$Property == "Silt",], aes(Laboratorybased, visNIRaided)) + theme_bw() + 
      labs(y = "Silt content, % (vis-NIR augmented)", x = "Silt content, % (Laboratory-based)") +
      geom_point (alpha = 0.5, size = 3) + geom_abline(intercept = 0, slope = 1, colour = "red") +
      xlim(0, 30) + ylim(0, 30) +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.y  = element_text(angle=90, vjust =0.5, hjust =0.5, size=14), legend.position="top", legend.title=NULL) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.x  = element_text(angle=0, vjust=0, size=14)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=18),
            legend.title = element_text(face="bold", colour = grey(0.2), size=18)) +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold",size = 16, colour = "black", angle = 0)) 

    d2silt + facet_wrap(~ Layer, nrow = 1)



    d2clay <- ggplot(pixpreds[pixpreds$Property == "Clay",], aes(Laboratorybased, visNIRaided)) + theme_bw() + 
      labs(y = "Clay content, % (vis-NIR augmented)", x = "Clay content, % (Laboratory-based)") +
      geom_point (alpha = 0.5, size = 3) + geom_abline(intercept = 0, slope = 1, colour = "red") +
      xlim(0, 70) + ylim(0, 70) +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.y  = element_text(angle=90, vjust =0.5, hjust =0.5, size=14), legend.position="top", legend.title=NULL) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.x  = element_text(angle=0, vjust=0, size=14)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=18),
            legend.title = element_text(face="bold", colour = grey(0.2), size=18)) +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold",size = 16, colour = "black", angle = 0)) 

    d2clay + facet_wrap(~ Layer, nrow = 1)

    d2ca <- ggplot(pixpreds[pixpreds$Property == "Ca",], aes(Laboratorybased, visNIRaided)) + theme_bw() + 
      labs(y = expression(paste("", Ca^"++", ", ", mmol[c], kg^-1, " (vis-NIR augmented)")), 
           x = expression(paste("", Ca^"++", ", ", mmol[c], kg^-1, " (Laboratory-based)"))) +
      geom_point (alpha = 0.5, size = 3) + geom_abline(intercept = 0, slope = 1, colour = "red") +
      xlim(0, 70) + ylim(0, 70) +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.y  = element_text(angle=90, vjust =0.5, hjust =0.5, size=14), legend.position="top", legend.title=NULL) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=18),
            axis.text.x  = element_text(angle=0, vjust=0, size=14)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=18),
            legend.title = element_text(face="bold", colour = grey(0.2), size=18)) +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold",size = 16, colour = "black", angle = 0)) 

    d2ca + facet_wrap(~ Layer, nrow = 1)

Define palette of color for the plot

    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

Clay maps

    clay_map <- ggplot(preditions_final[preditions_final$method != "Differences between maps",], 
                       aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Clay)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.y  = element_text(angle=0, vjust =0.5, hjust =0.5, size=16), 
            legend.title = element_text(face = "bold", colour = "black", size=25)) +
      guides(fill=guide_colourbar(title="Clay content, % ")) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_blank()) 

    clay_map_diff <- ggplot(preditions_final[preditions_final$method == "Differences between maps",], 
                       aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Clay)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(), 
            legend.title = element_text(face = "bold", colour = "black", size = 25)) +
      guides(fill=guide_colourbar(title="Difference, % ")) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_text(face = "bold", size = 25, colour = "black", angle = -90)) 
    clay_map_diff

Save the plots in pdf format

    pdf(file = "Map_Clay.pdf", width = 20, height = 11)
    grid.draw(cbind(ggplotGrob(clay_map), 
                    ggplotGrob(clay_map_diff + theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))), 
                    size = "last"),
              recording=F)
    dev.off()

Silt maps

    silt_map <- ggplot(preditions_final[preditions_final$method != "Differences between maps",], 
                       aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Silt)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.y  = element_text(angle=0, vjust =0.5, hjust =0.5, size=16), 
            legend.title = element_text(face = "bold", colour = "black", size=25)) +
      guides(fill=guide_colourbar(title="Silt content, % ")) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_blank()) 

    silt_map_diff <- ggplot(preditions_final[preditions_final$method == "Differences between maps",], 
                            aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Silt)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(), 
            legend.title = element_text(face = "bold", colour = "black", size = 25)) +
      guides(fill=guide_colourbar(title="Difference, % ")) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_text(face = "bold", size = 25, colour = "black", angle = -90)) 
    silt_map_diff

Save the plots in pdf format

    pdf(file = "Map_Silt.pdf", width = 20, height = 11)
    grid.draw(cbind(ggplotGrob(silt_map), 
                    ggplotGrob(silt_map_diff + theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))), 
                    size = "last"),
              recording=F)
    dev.off()

Sand maps

    sand_map <- ggplot(preditions_final[preditions_final$method != "Differences between maps",], 
                       aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Sand)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.y  = element_text(angle=0, vjust =0.5, hjust =0.5, size=16), 
            legend.title = element_text(face = "bold", colour = "black", size=25)) +
      guides(fill=guide_colourbar(title="Sand content, % ")) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_blank()) 

    sand_map_diff <- ggplot(preditions_final[preditions_final$method == "Differences between maps",], 
                            aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Sand)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(), 
            legend.title = element_text(face = "bold", colour = "black", size = 25)) +
      guides(fill=guide_colourbar(title="Difference, % ")) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_text(face = "bold", size = 25, colour = "black", angle = -90)) 
    sand_map_diff

Save the plots in pdf format

    pdf(file = "Map_Sand.pdf", width = 20, height = 11)
    grid.draw(cbind(ggplotGrob(sand_map), 
                    ggplotGrob(sand_map_diff + theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))), 
                    size = "last"),
              recording=F)
    dev.off()

    pdf(file = "Map_Sand.pdf", width = 20, height = 11)
    grid.draw(rbind(cbind(ggplotGrob(sand_map), 
                          ggplotGrob(sand_map_diff + theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))), 
                          size = "last")))
                   
    dev.off()

Ca maps

    ca_map <- ggplot(preditions_final[preditions_final$method != "Differences between maps",], 
                     aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Ca)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.y  = element_text(angle=0, vjust =0.5, hjust =0.5, size=16), 
            legend.title = element_text(face = "bold", colour = "black", size=25)) +
      guides(fill=guide_colourbar(title=expression(paste("", Ca^"++", ", ", mmol[c], kg^-1, " ")))) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_blank()) 


    ca_map_diff <- ggplot(preditions_final[preditions_final$method == "Differences between maps",], 
                          aes(POINT_X, POINT_Y)) + geom_tile(aes(fill = Ca)) + 
      facet_grid(layer ~ method) +
      scale_fill_gradientn(colours = myPalette(30)) + coord_equal() + theme_bw() + 
      theme(legend.position = "top") +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(), 
            legend.title = element_text(face = "bold", colour = "black", size = 25)) +
      guides(fill=guide_colourbar(title=expression(paste("Difference", ", ", mmol[c], kg^-1, " ")))) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.2), size=25),
            axis.text.x  = element_text(angle = 0, vjust=0, size=16)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.2), size=15),
            legend.key.size = unit(0.95, "cm")) +
      labs(y = "Northings (m)", x = "Eastings (m)") +
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold", size = 25, colour = "black", angle = 0),
            strip.text.y = element_text(face = "bold", size = 25, colour = "black", angle = -90)) 
    ca_map_diff

Save the plot in pdf format

    pdf(file = "Map_Ca.pdf", width = 20, height = 11)
    grid.draw(cbind(ggplotGrob(ca_map), 
                    ggplotGrob(ca_map_diff + theme(plot.margin = unit(c(0,0.5,0,0.5), "cm"))), 
                    size = "last"),
              recording=F)
    dev.off()

Save the results

    save.image("spatial_modeling.RData")
