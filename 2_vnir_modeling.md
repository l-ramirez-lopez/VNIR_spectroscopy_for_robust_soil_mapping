========================================================

Refers to section vis-NIR modelling and predictions (modelling)

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
    if (!require(rgdal))
      install.packages("rgdal") # version 1.2-7
    if (!require(sp))
      install.packages("sp") # version 1.2-4
    if (!require(ggplot2))
      install.packages("ggplot2") # version 2.2.1
    if (!require(maptools))
      install.packages("maptools") # version 0.9-2
    if (!require(plyr))
      install.packages("plyr") # version 1.8.4
    if (!require(wesanderson))
      install.packages("wesanderson")
    if (!require(beepr))
      install.packages("beepr") # version 1.2
    options(warn=0)

Define the working directory

    workingd <- "your/working/directory/"
    setwd(workingd)

Organize the data
=================

    data <- read.table("SoilNIRSaoPaulo.txt", header = TRUE, check.names = FALSE, sep ="\t")

Minimum wavelength in the dataset

    mwav <- 502
    spc <- log(1/data[,which(colnames(data) == mwav):ncol(data)])
    data <- data[,-c(which(colnames(data) == mwav):ncol(data))]
    data$spc <- spc

Apply standard normal variate

    data$spc_snv <- standardNormalVariate(data$spc)
    rm(spc)

    nr <- range(as.numeric(colnames(data$spc)))

Split the dataset

    data$alr_silt <- data$Silt/data$Sand
    data$alr_clay <- data$Clay/data$Sand

    valida <- data[data$set == "validation",]
    data <- data[data$set == "cal_candidate",]

Compress the data

    pcaall <- orthoProjection(Xr = data$spc_snv, X2 = NULL, 
                              Yr = NULL, 
                              method = "pca", pcSelection = list("var", 0.01), 
                              center = TRUE, scaled = FALSE, cores = 1)

Standardize the scores

    pcaall$scores.std <- sweep(pcaall$scores, MARGIN = 2, STATS = pcaall$sc.sdv, FUN = "/")

Select the final calibration samples
====================================

Repeat the sampling several times and select the one with the minimum
mds (ramirez-lopez et al., 2014)

Optimal calibration sample set size

    ocss <- 180

Compute the maximum and minimum of each score for the limits of the
density estimations

    max.sc <- colMins(pcaall$scores.std)
    min.sc <- colMaxs(pcaall$scores.std)

Compute the mean and standard deviation of each score for the
comparisons with the samples

    mean.sc <- colMeans(pcaall$scores.std)
    sd.sc <- colSds(pcaall$scores.std)

Number of points in the density distribution

    nxdens <- 2000

Matrix where the density values will be stored

    i= 1
    sc.dens <- data.frame(seq(min.sc[i], max.sc[i], length = nxdens), 
                          matrix(NA, nxdens, length(min.sc)))
    colnames(sc.dens) <- c("x", paste("densc", 1:length(min.sc), sep = ""))

Estimate the density distribution of each component

    d.bandwidths <- rep(NA, length(min.sc))
    for(i in 1:length(min.sc)){
      i.pdf <- density(pcaall$scores.std[,i], 
                       bw = "nrd0", 
                       n = nxdens, from = min.sc[i], to = max.sc[i], 
                       kernel = "gaussian")
      d.bandwidths[[i]] <- i.pdf$bw
      sc.dens[,i+1] <- i.pdf$y
    }


    reps <- 10
    results.clhs.ocss <- rep(NA, length = reps)
    calidx <- NULL
    for(i in 1:reps)
    {
      set.seed(round(i*exp(4)))
      i.calidx <- clhs(as.data.frame(pcaall$scores.std),
                       size = ocss, iter = 10000)

      j.sc.dens <- msd.sc  <- sc.dens
      for(j in 1:length(min.sc)){
        j.sc.dens[,j+1] <- density(pcaall$scores.std[i.calidx,j],
                                   bw = d.bandwidths[j],
                                   n = nxdens, from = min.sc[j], to = max.sc[j],
                                   kernel = "gaussian")$y
      }
      results.clhs.ocss[i] <- mean(colMeans((j.sc.dens[,-1] - sc.dens[,-1])^2, na.rm = T))

      calidx[[i]] <- i.calidx
    }

    plot(results.clhs.ocss)

Best group

    calibration.idx <- calidx[[6]]
    cal_smpls <- as.character(data$Id[calibration.idx])

Save the results

    save.image("vnir_modeling.RData")

Divide the data in calibration and validation sets
==================================================

    train <- data[as.character(data$Id) %in%  cal_smpls, ]
    pred <- data[!(as.character(data$Id) %in%  c(cal_smpls)), ]

    train$layer <- as.factor(substr(train$Id, 1, 1))
    valida$layer <- as.factor(substr(valida$Id, 1, 1))
    pred$layer <- as.factor(substr(pred$Id, 1, 1))

    train$Id <- as.factor(as.character(train$Id))
    valida$Id <- as.factor(as.character(valida$Id))
    pred$Id <- as.factor(as.character(pred$Id))

Additive log ratio transformation using Sand as the denominator
retreives the best models

Use the spectrum-based learner (sbl) algorithm:

Ditance computations: principal components distance

Usage of the distance information: source of additional predictors

Regression method: gaussian process

    ctrl <- mblControl(sm = "pls", pcSelection = list("opc", 40), 
                       valMethod = c("NNv"), 
                       returnDiss = TRUE,
                       scaled = FALSE, center = TRUE)

Define the threshold distances to test Note: the number of threshold
distances tested here is large, it could be reduced

    diss2test <- seq(0.3, 1.5, by = 0.05)

Define the minimum and maximum number of neighbors to be included in the
model.

Example: If with a given threshold distance only 70 neighbors are
selected, then the function will be fored to take the especified minimum
number of neighbors to be included in the model. In this case 120. In
this case we also set that the maximum number of neighbors that can be
included can be all the samples in the calibration set.

    kminmax <- c(120, nrow(train$spc))

Optimize the number threshold distance for neighbor selection and
predict Ca for valication

    sbl_ca <- mbl(Yr = train$Ca, 
                  Xr = train$spc,
                  Xu = valida$spc,
                  mblCtrl = ctrl,
                  group = substr(train$Id, 2, 100),
                  dissUsage = "none",
                  k.diss = diss2test, 
                  k.range = kminmax,
                  pls.c = c(5,20),
                  method = "wapls1")

    idx.best.ca <- which.min(sbl_ca$nnValStats$st.rmse)
    best.kdiss.ca <- sbl_ca$nnValStats$k.diss[idx.best.ca]

    pred_valCa <- getPredictions(sbl_ca)[,idx.best.ca]

    valR2Ca <- cor(valida$Ca, pred_valCa)^2
    valrmseCa <- (mean((valida$Ca - pred_valCa)^2))^0.5
    valmeCa <- mean(valida$Ca - pred_valCa)

    rownames(valida) <- seq(1,nrow(valida), by=1)
    var.res.Ca_a <- var(valida[valida$layer=='A','Ca'] - pred_valCa[as.numeric(rownames(valida[valida$layer=='A',]))])
    var.res.Ca_b <- var(valida[valida$layer=='B','Ca'] - pred_valCa[as.numeric(rownames(valida[valida$layer=='B',]))])
    var.res.Ca <- var(valida$Ca - pred_valCa)

    valR2Ca;valrmseCa;valmeCa

Number of samples that are predicted by local regression (instead of
using all the neighbors, i.e. global)

    sum(sbl_ca$results[[idx.best.ca]]$k.org != nrow(train$spc))

Optimize the number threshold distance for neighbor selection and
predict alr silt for validation

Note: the number of threshold distances tested here is large, it could
be reduced

    sbl_alrsilt <- mbl(Yr = log(train$Silt/train$Sand), 
                       Xr = train$spc, 
                       Xu = valida$spc,
                       mblCtrl = ctrl,
                       group = substr(train$Id, 2, 100),
                       dissUsage = "none",
                       k.diss = diss2test, 
                       k.range = kminmax,
                       pls.c = c(5,20),
                       method = "wapls1")

    idx.best.alrsilt <- which.min(sbl_alrsilt$nnValStats$st.rmse)
    best.kdiss.alrsilt <- sbl_alrsilt$nnValStats$k.diss[idx.best.alrsilt]

    pred_val_alrsilt <- getPredictions(sbl_alrsilt)[,idx.best.alrsilt]

    rownames(valida) <- seq(1,nrow(valida), by=1)

    var.res.alr_silt_a <- var(valida[valida$layer=='A','alr_silt'] - pred_val_alrsilt[as.numeric(rownames(valida[valida$layer=='A',]))])
    var.res.alr_silt_b <- var(valida[valida$layer=='B','alr_silt'] - pred_val_alrsilt[as.numeric(rownames(valida[valida$layer=='B',]))])
    var.res.alrsilt <- var(valida$alr_silt - pred_val_alrsilt)

Number of samples that are predicted by local regression (instead of
using all the neighbors, i.e. global)

    sum(sbl_alrsilt$results[[idx.best.alrsilt ]]$k.org != nrow(train$spc))

Optimize the number threshold distance for neighbor selection and
predict alr silt for validation Note: the number of threshold distances
tested here is large,it could be reduced

    sbl_alrclay <- mbl(Yr = log(train$Clay/train$Sand), 
                       Xr = train$spc, 
                       Xu = valida$spc,
                       mblCtrl = ctrl,
                       group = substr(train$Id, 2, 100),
                       dissUsage = "none",
                       k.diss = diss2test, 
                       k.range = kminmax,
                       pls.c = c(5,20),
                       method = "wapls1")

    idx.best.alrclay <- which.min(sbl_alrclay$nnValStats$st.rmse)
    best.kdiss.alrclay <- sbl_alrclay$nnValStats$k.diss[idx.best.alrclay]

    pred_val_alrclay <- getPredictions(sbl_alrclay)[,idx.best.alrclay]


    rownames(valida) <- seq(1,nrow(valida), by=1)

    var.res.alr_clay_a <- var(valida[valida$layer=='A','alr_clay'] - pred_val_alrclay[as.numeric(rownames(valida[valida$layer=='A',]))])
    var.res.alr_clay_b <- var(valida[valida$layer=='B','alr_clay'] - pred_val_alrclay[as.numeric(rownames(valida[valida$layer=='B',]))])
    var.res.alrclay <- var(valida$alr_clay - pred_val_alrclay)

Number of samples that are predicted by local regression (instead of
using all the neighbors, i.e. global)

    sum(sbl_alrclay$results[[idx.best.alrclay]]$k.org != nrow(train$spc))

    valSilt_alr <- 100 * exp(pred_val_alrsilt)/(1+exp(pred_val_alrclay)+exp(pred_val_alrsilt))
    valClay_alr <- 100 * exp(pred_val_alrclay)/(1+exp(pred_val_alrclay)+exp(pred_val_alrsilt))
    valSand_alr <- 100 * 1/(1+exp(pred_val_alrclay)+exp(pred_val_alrsilt))


    valR2Silt_alr <- (cor(valida$Silt, valSilt_alr))^2
    valR2Sand_alr <- (cor(valida$Sand, valSand_alr))^2
    valR2Clay_alr <- (cor(valida$Clay, valClay_alr))^2

    valrmseSilt_alr <- mean((valida$Silt - valSilt_alr)^2)^0.5
    valrmseSand_alr <- mean((valida$Sand - valSand_alr)^2)^0.5
    valrmseClay_alr <- mean((valida$Clay - valClay_alr)^2)^0.5

    valmeSilt_alr <- mean(valida$Silt - valSilt_alr)
    valmeSand_alr <- mean(valida$Sand - valSand_alr)
    valmeClay_alr <- mean(valida$Clay - valClay_alr)

Predict Ca in the prediction set based on the optimized threshold
distance for neighbor selection

    sbl_ca_pred <- mbl(Yr = train$Ca, 
                       Xr = train$spc, 
                       Xu = pred$spc,
                       mblCtrl = ctrl,
                       group = substr(train$Id, 2, 100),
                       dissUsage = "none",
                       k.diss = best.kdiss.ca, 
                       k.range = c(120, 180),
                       pls.c = c(5,20),
                       method = "wapls1")

    pred_predCa <- getPredictions(sbl_ca_pred)[,1]

Number of samples that are predicted by local regression (instead of
using all the neighbors, i.e. global)

    sum(sbl_ca_pred$results[[1]]$k.org != nrow(train$spc))

Predict alr Silt in the prediction set based on the optimized threshold
distance for neighbor selection

    sbl_alrsilt_pred <- mbl(Yr = log(train$Silt/train$Sand), 
                            Xr = train$spc, 
                            Xu = pred$spc,
                            mblCtrl = ctrl,
                            group = substr(train$Id, 2, 100),
                            dissUsage = "none",
                            k.diss = best.kdiss.alrsilt, 
                            k.range = c(120, 180),
                            pls.c = c(5,20),
                            method = "wapls1")

    pred_pred_alrsilt <- getPredictions(sbl_alrsilt_pred)[,1]

Number of samples that are predicted by local regression (instead of
using all the neighbors, i.e. global)

    sum(sbl_alrsilt_pred$results[[1]]$k.org != nrow(train$spc))

Predict alr Clay in the prediction set based on the optimized threshold
distance for neighbor selection

    sbl_alrclay_pred <- mbl(Yr = log(train$Clay/train$Sand), 
                            Xr = train$spc, 
                            Xu = pred$spc,
                            mblCtrl = ctrl,
                            group = substr(train$Id, 2, 100),
                            dissUsage = "none",
                            k.diss = best.kdiss.alrclay, 
                            k.range = c(120, 180),
                            pls.c = c(5,20),
                            method = "wapls1")

    pred_pred_alrclay <- getPredictions(sbl_alrclay_pred)[,1]

Number of samples that are predicted by local regression (instead of
using all the neighbors, i.e. global)

    sum(sbl_alrclay_pred$results[[1]]$k.org != nrow(train$spc))

Predict

    predSilt_alr <- 100 * exp(pred_pred_alrsilt)/(1+exp(pred_pred_alrclay)+
                                                    exp(pred_pred_alrsilt))
    predClay_alr <- 100 * exp(pred_pred_alrclay)/(1+exp(pred_pred_alrclay)+
                                                    exp(pred_pred_alrsilt))
    predSand_alr <- 100 * 1/(1+exp(pred_pred_alrclay)+exp(pred_pred_alrsilt))

Store the results

    results_spec_modeling <- data.frame(Property = c("Sand", "Silt", "Clay", "Ca"), 
                                        R2 = rep(NA, 4), 
                                        RMSE = rep(NA, 4),
                                        ME = rep(NA, 4))
                                
    results_spec_modeling$R2 <- c(valR2Sand_alr, valR2Silt_alr, 
                                  valR2Clay_alr, valR2Ca)
    results_spec_modeling$RMSE <- c(valrmseSand_alr, valrmseSilt_alr, 
                                    valrmseClay_alr, valrmseCa)
    results_spec_modeling$ME <- c(valmeSand_alr, valmeSilt_alr, 
                                    valmeClay_alr, valmeCa)

Write table of the results

    write.table(results_spec_modeling, file = "results_spec_modeling.txt", 
                row.names = FALSE, sep = "\t")

Plot the results

    locsamples <- rbind(data.frame(ns = 
                                as.numeric(as.character(sbl_alrsilt_pred$results[[1]]$k)), 
                                   Property = "alrsilt"),
                        data.frame(ns = 
                                as.numeric(as.character(sbl_alrclay_pred$results[[1]]$k)), 
                                   Property = "alrclay"),
                        data.frame(ns = 
                                as.numeric(as.character(sbl_ca_pred$results[[1]]$k)), 
                                   Property = "Ca"))


    sum(sbl_ca_pred$results[[1]]$k.org == nrow(train$spc))/nrow(pred)
    sum(sbl_alrsilt_pred$results[[1]]$k.org == nrow(train$spc))/nrow(pred)
    sum(sbl_alrclay_pred$results[[1]]$k.org == nrow(train$spc))/nrow(pred)

    ggplot(locsamples, aes(x = ns, fill = Property))+
      geom_histogram(aes(y=..density..),
                     alpha=0.5,position='identity',binwidth=10)

Prepare the vis-NIR-aided database for spatial modelling
========================================================

    specaidaddb <- rbind(data.frame(train[,c("Id", "POINT_X", "POINT_Y", 
                                             "layer", "Clay", "Sand", "Silt", "Ca")], 
                                    set = "train", 
                                    specClay_alr = train$Clay, 
                                    specSand_alr = train$Sand, 
                                    specSilt_alr = train$Silt,
                                    specCa = train$Ca), 
                         
                         data.frame(valida[,c("Id", "POINT_X", "POINT_Y", 
                                              "layer", "Clay", "Sand", "Silt", "Ca")], 
                                    set = "validation", 
                                    specClay_alr = valida$Clay, 
                                    specSand_alr = valida$Sand, 
                                    specSilt_alr = valida$Silt,
                                    specCa = valida$Ca),  
                         
                         data.frame(pred[,c("Id", "POINT_X", "POINT_Y", 
                                            "layer", "Clay", "Sand", "Silt", "Ca")], 
                                    set = "prediction", 
                                    specClay_alr = predClay_alr, 
                                    specSand_alr = predSand_alr, 
                                    specSilt_alr = predSilt_alr,
                                    specCa = pred_predCa))

Write table of the vis-NIR database

    write.table(specaidaddb, file = "spec_aideddb.txt", row.names = F, sep = "\t")

Plot and save the data

    wav <- as.numeric(colnames(train$spc))
    spca <- 1/exp(data$spc[substr(data$Id, 1, 1) == "A", ])
    spcb <- 1/exp(data$spc[substr(data$Id, 1, 1) == "B", ])

    reflectance.all <- rbind(data.frame(Reflectance = 100 * as.vector(t(spca)),
                                        "Wavelength (nm)" = rep(wav, nrow(spca)),
                                        Layer = "Layer A (0 - 0.2 m)",
                                        Sample = paste("Sample_", 
                                        rep(1:nrow(spca), each = ncol(spca)), sep = "")),
                             data.frame(Reflectance = 100 * as.vector(t(spcb)),
                                        "Wavelength (nm)" = rep(wav, nrow(spcb)),
                                        Layer = "Layer B (0.8 - 1.0 m)",
                                        Sample = paste("Sample_", 
                                        rep(1:nrow(spcb), each = ncol(spcb)), sep = "")))

    p.tmp.all <- ggplot(reflectance.all, 
                        aes(x = Wavelength..nm., y = Reflectance, shape = Sample)) + 
      theme_bw() +
      geom_line(size = 0.5, colour =  rgb(0.1,0.1,0.1,0.20), linetype = 1) +
      theme(axis.title.y = element_text(face="bold", colour = grey(0.25), size = 60),
            axis.text.y  = element_text(angle=90, vjust =0.5, hjust =0.5, size = 40), 
            legend.position="top", legend.title=NULL) +
      theme(axis.title.x = element_text(face="bold", colour = grey(0.25), size = 60),
            axis.text.x  = element_text(angle=0, vjust=0, size = 40)) +
      theme(legend.text = element_text(face="bold", colour = grey(0.25), size = 1),
            legend.title = element_text(face="bold", colour = grey(0.25), size = 1)) + 
      theme(strip.background = element_rect(fill = "grey"), 
            strip.text.x = element_text(face = "bold",size = 60, colour = "black", 
                                        angle = 0)) +
      labs(y = "Reflectance, %", x = "Wavelength (nm)") +
      ylim(0, 50)

    p.tmp.all + facet_grid(.~Layer) + theme(plot.margin = unit(c(1,1,1,1), "cm"))


    pdf(file = "spectra_layers.pdf", width = 30, height = 15)
    p.tmp.all + facet_grid(.~Layer)  + theme(plot.margin = unit(c(1,1,1,1), "cm"))
    dev.off()

    cormat <- (cor(cbind(data[,c("Sand", "Silt", "Clay", "Ca")], 
                         alr_silt = data$Silt/data$Sand, 
                         alr_Clay = data$Clay/data$Sand, 
                         mean_reflectance = rowMeans(data$spc))))

    write.table(round(cormat, 2), file = "correlation_matrix.txt", row.names = F, sep = "\t")

Save the variance of the spectroscopic model error

    save(var.res.Ca, var.res.Ca_a, var.res.Ca_b, var.res.alr_clay_a, var.res.alr_clay_b, var.res.alrclay,var.res.alr_silt_a, var.res.alr_silt_b, var.res.alrsilt, file = "var_spec_error.Rdata")
