# Robust soil mapping at the farm scale with vis-NIR spectroscopy 
_Leo Ramirez-Lopez & Alex Wadoux_

<hr>
This repository is dedicated to our paper entitled ['Robust soil mapping at the farm scale with vis-NIR spectroscopy'](http://onlinelibrary.wiley.com/doi/abs/10.1111/ejss.12752) published in European Journal of Soil Science. 
In the spirit of reproducible research, here you will find all the the data and computational code we used for conducting the research associated to this paper. If it catches your interest, you can go through the code presented here (mostly using copy and paste &#x263A;) and reproduce all the results reported in the paper. If you find that the data, parts of the code or some of the ideas here are useful for your own work, please cite any of these as follows:


<ul><li style="color: #FFFFFF;">[Ramirez-Lopez, L., Wadoux, A.M.C., Franceschini, M. H. D., Terra, F.S., Marques, K.P.P., Sayão, V.M., Demattê, J.A.M. 2018.Robust soil mapping at the farm scale with vis-NIR spectroscopy. European Journal of Soil Science. doi: 10.1111/ejss.12752.](http://onlinelibrary.wiley.com/doi/abs/10.1111/ejss.12752)</li></ul>


If you want to check the code and try to reproduce what we present in the paper, we suggest to use [R >= 3.4.4](https://cran.r-project.org/) and also use [RStudio (>= 1.1.442)](https://www.rstudio.com/products/RStudio) which is a really good integrated development environment (IDE) for R. 

# The Barra Bonita dataset
<p>__License__: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)</p>
This is a soil dataset that comes from a soil survey conducted in a farm located in the municipality of [Barra Bonita](https://en.wikipedia.org/wiki/Barra_Bonita,_S%C3%A3o_Paulo) in the state of São Paulo (Brazil). This dataset comprises 910 soil samples collected at 458 locations (at two depths). The area covered by the survey is 473 ha.  
The soil dataset used in the study is provided [here](https://github.com/l-ramirez-lopez/VNIR_spectroscopy_for_robust_soil_mapping/raw/master/BarraBonita.txt) in a table format which comprises the following headers (variable names):
<ul>
<li>__`Nr`__: An arbitrary number.</li>
<li>__`ID`__: The sample identifier which is given by a letter followed by a number. The letter designates the depth at which the soil the soil sample was collected. The letter A stands for the samples collected at a depth of 0-0 .2 m (458 samples) and letter B stands for samples collected at a depth of 0.8-1.0 m (452 samples). The number designates an arbitrary number assigned to the sampling location (same profile).</li>
<li>__`POINT_X`__ and	__`POINT_Y`__: These two variables are the geographical coordiates of the samples.</li>
<li>__`Sand`__, __`Silt`__ and __`Clay`__:	These variables are the sand, silt and sand contnets (soil particle-size fractions) of the samples (in percetnage).</li>
<li>__`Ca`__: This variable is the exchangeable calcium content (in mmol~c~ kg^-1^).</li>
<li>__`set`__: Is a variable with two categories (`validation` and `cal_candidate`). The category `validation` indicate a subset samples (227) that were selected by random selection of 115 sampling locations. These samples were used in our work to conduct the validation of the different models. 
<li>__`502-2338`__: These are 307 spectral reflectance variables corresponding to the visible (vis) and near-infrared (NIR) region of the electromagnetic spectrum. Each variable name is a number which represents the wavelength (in nm). This data has been resampled from its original resolution to a new resolution of 6 nm.</li>
</ul>
Additional details of this dataset are provided in our paper.

# The code
The code presented here is divided in three main parts according to the sections in the paper as follows:
<ul>
<li> Calibration sampling: provides the code used in the section 'Calibration samples'.</li>
<li> Modeling of the vis-NIR spectra: provides the code used in the sections 'Transformation of the particle-size data' and  'Vis-NIR modelling and predictions'. </li>
<li> Spatial modeling: provides the codes used in the section 'Spatial modelling of the vis-NIR augmented data' of the paper.</li>
</ul>

Before you start, make sure you get all the necessary packages installed and download the Barra Bonita dataset.

Install the required packages

```r
## (optional) set the language to English
Sys.setenv(language = "EN")

## List of the R packages that you will need to have installed
requiredp <- c("resemble",     # version 1.2.2
               "prospectr",    # version 0.1.3
               "clhs",         # version 0.5-6
               "matrixStats",  # version 0.51.0
               "doParallel",   # version 1.0.10
               "ggplot2",      # version 2.2.1
               "rgdal",        # version 1.2-7
               "sp",           # version 1.2-4
               "maptools",     # version 0.9-2
               "plyr",         # version 1.8.4
               "raster",       # version 2.5-8
               "gridExtra",    # version 2.2.1
               "RColorBrewer", # version 1.1-2
               "grid",         # version 3.3.2
               "cowplot",      # version 1.1-2
               "extrafont"     # version 3.3.2
               )
## Check which packages are not yet installed in your R
ptoinstall <- requiredp[!requiredp %in% installedp]

## If there are are packages not yet installed .. then install them 
if(length(ptoinstall) > 0)
  install.packages(ptoinstall)

## Load all the necessary packages
lapply(requiredp, require, character.only = TRUE)
```

Define the working directory

```r
workingd <- "C:/your/working/directory"
setwd(workingd)
```
Download the Barra Bonita dataset
```r
## File URL
dataurl <- "https://github.com/l-ramirez-lopez/VNIR_spectroscopy_for_robust_soil_mapping/raw/master/BarraBonita.txt"

## download the file and store it as "BarraBonita.txt" in your working directory
download.file(dataaurl, destfile = "BarraBonita.txt", method = "auto")
```



## 1. Selecting the candidate samples for calibration and the validation set

Read the data

```r
data <- read.table("BarraBonita.txt", header = TRUE, check.names = FALSE, sep = "\t")
```

Organize the data
```r
## Minimum wavelength in the dataset
mwav <- 502

## Select the spectral variables and transform them to apparent absorbance
## and put the spectral data into a separate object
spc <- log(1/data[,which(colnames(data) == mwav):ncol(data)])

## Overwrite the data with the non-spectral variables only 
data <- data[,-c(which(colnames(data) == mwav):ncol(data))]

## put the spectral data back into the data
data$spc <- spc

## remove the unnecessary object
rm(spc)
```

Apply standard normal variate

```r
data$spc_snv <- standardNormalVariate(data$spc)

```

Remove the validation samples from the analysis

```r
valida <- data[data$set == "validation",]
data <- data[data$set == "cal_candidate",]
```

__1.1 Perform a principal component analysis__
<p>Project the spectra onto the principal components (PCs) space.</p>

```r
pcaall <- orthoProjection(Xr = data$spc_snv, X2 = NULL, 
                          Yr = NULL, 
                          method = "pca", pcSelection = list("var", 0.01), 
                          center = TRUE, scaled = FALSE, cores = 1)

## Standardize the scores
pcaall$scores.std <- sweep(pcaall$scores, MARGIN = 2, STATS = pcaall$sc.sdv, FUN = "/")
```
__1.2 Compute the probability density functions (pdfs) of each PC__

<p>The probability density functions (pdfs) estimated here will be used as references to compare the pdfs computed with the different subsets of different sizes.</p>


```r
## First, compute the maximum and minimum of each score for the limits of the density estimations
max.sc <- colMins(pcaall$scores.std)
min.sc <- colMaxs(pcaall$scores.std)

## Compute the mean and standard deviation of each score for the comparisons with the samples
mean.sc <- colMeans(pcaall$scores.std)
sd.sc <- colSds(pcaall$scores.std)

## Indicate a number of poits to get the density function
nxdens <- 500

## Define the matrix where the density values will be stored
ix <- 1
sc.dens <- data.frame(seq(min.sc[ix], max.sc[ix], length = nxdens), 
                      matrix(NA, nxdens, length(min.sc)))
colnames(sc.dens) <- c("x", paste("densc", 1:length(min.sc), sep = ""))

## Finally, Estimate the density function of each component
d.bandwidths <- rep(NA, length(min.sc))
names(d.bandwidths) <- colnames(pcaall$scores.std)
for(i in 1:length(min.sc)){
  idsty <- density(pcaall$scores.std[,i], 
                   bw = "nrd0", 
                   n = nxdens, from = min.sc[i], to = max.sc[i], 
                   kernel = "gaussian")
  sc.dens[,i+1] <- idsty$y
  d.bandwidths[i] <- idsty$bw
}
```

__1.3 Select subsets of different sizes and compare their pdfs against the reference pdfs__

The following steps might take a while...
```r
## Create a vector containing the different sample sizes to test
css <- seq(10, 400, by = 10)

## Create a matrix to store the results of the comparisons
results.clhs <-  data.frame(css = css, 
                            msd = rep(NA, length(css)),
                            mndiff = rep(NA, length(css)),
                            sddiff = rep(NA, length(css)))


## Indicate the number of repetitions to be used in order to get 
## reliable estimates of the msd as a function of the sample set size
sreps <- 10

for(k in 1:sreps){
  
  cat("--------------- Repetition", k, "---------------\n")
  
  results.clhs[,-1] <- NA
  
  fn <- paste("5pcs_resultsclhs_rep", k,".txt", sep = "")
  
  if(fn %in% list.files()){
    results.clhs <- read.table(fn, header = T, sep = "\t")
  }
  
  iter.p <- 1 + sum(rowSums(!is.na(results.clhs)) == ncol(results.clhs))
  
  for(i in iter.p:length(css)){
    
    cat("------ Sample size:", css[i], "; Repetition", k, "------\n")
    
    set.seed(k)
    
    ## Sample with conditioned latin hypercube sampling algorithm
    i.calidx <- clhs(as.data.frame(pcaall$scores.std),
                     size = css[i], iter = 10000)
    
    j.sc.dens <- msd.sc  <- sc.dens 
    for(j in 1:length(min.sc)){
      
      ## Compute the pdf of the jth PC for the ith sample size 
      ## use the same bandwidth (bw) as in the whole set of candidates
      j.sc.dens[,j+1] <- density(pcaall$scores.std[i.calidx,j], 
                                 bw = d.bandwidths[j], 
                                 n = nxdens, from = min.sc[j], to = max.sc[j], 
                                 kernel = "gaussian")$y
    }
    ## get the msd
    results.clhs$msd[i] <- mean(colMeans((j.sc.dens[,-1] - sc.dens[,-1])^2, na.rm = T))
    
    ## (optional) compare the means of the sample against the 
    ## mean of the whole set
    results.clhs$mndiff[i] <- mean(abs(colMeans(pcaall$scores.std[i.calidx,])))
    
    ## (optional) compare the standard deviations of the sample 
    ## against the mean of the whole set
    results.clhs$sddiff[i] <- mean(abs(colSds(pcaall$scores.std[i.calidx,]) - 1))
    
    
    ## write the results into a file 
    write.table(results.clhs, 
                file = fn,
                row.names = FALSE, sep = "\t")
    
    print(results.clhs[1:i,])
    cat("\n")
  }
}

## Read all the results from the iterations and compute the average 
## of the msds
nmsrepsclhs <- paste("5pcs_resultsclhs_rep", 1:10, ".txt", sep = "")
final.clhs <- 0
for(i in nmsrepsclhs){
  iter <- which(i == nmsrepsclhs)
  results.clhs <- read.table(i, header = T, sep = "\t")
  results.clhs$mndiff <- abs(results.clhs$mndiff)
  final.clhs <- final.clhs + results.clhs
  if(iter == length(nmsrepsclhs))
  {
    final.clhs <- final.clhs/iter
    write.table(final.clhs, file = "6pcs_final.clhs.txt", sep = "\t", row.names = FALSE)
  }
}

## Read all the results from the iterations and compute the standard 
## deviation of the msds
final.clhs_sd <- 0
for(i in nmsrepsclhs){
  iter <- which(i == nmsrepsclhs)
  results.clhs <- read.table(i, header = T, sep = "\t")
  results.clhs$mndiff <- abs(results.clhs$mndiff)
  final.clhs_sd <- (results.clhs - final.clhs_sd)^2
  if(iter == length(nmsrepsclhs))
  {
    final.clhs_sd <- (final.clhs_sd/iter)^0.5
  }
}

final.clhs_plot <- data.frame(final.clhs[,1:2], 
                              sd_lower = final.clhs[,2] - final.clhs_sd[,2],
                              sd_upper = final.clhs[,2] + final.clhs_sd[,2])

## Plot the results using ggplot
p.tmp <- ggplot(final.clhs_plot) + geom_line(aes(x = css, msd, colour = "msd")) + 
  theme_bw() + 
  theme(axis.title.y = element_text(face= "bold.italic", colour = grey(0.2), size=18),
        axis.text.y  = element_text(angle=0, vjust =0.5, hjust =0.5, size=14), 
        legend.title = element_text(colour = "white", size=20)) +
  theme(axis.title.x = element_text(face= "bold", colour = grey(0.2), size=18),
        axis.text.x  = element_text(angle = 0, vjust=0, size=14)) +
  theme(legend.position = "none") + 
  theme(legend.text = element_text(face="bold", colour = grey(0.2), size=18)) +
  labs(y = "msd", x = "Calibration set size") +
  theme(strip.background = element_rect(fill = "grey"), 
        strip.text.x = element_text(size = 16, colour = "black", angle = 0)) 

## Save the plot
pdf(file = "calibration_set_size.pdf", width = 8, height = 6)
p.tmp + geom_ribbon(aes(ymin=sd_lower, ymax=sd_upper, x = css, colour = "bands"), alpha = 0.2) +
  scale_colour_manual(name = '', values = c("bands" = NA, "msd" = "black"))
dev.off()

final.clhs <- read.table("6pcs_final.clhs.txt", header = TRUE, check.names = FALSE, sep ="\t")
```

The plot that was generated above shows the mean squared Euclidean distance (msd) values corresponding to the comparisons between the estimates of probability density functions (pdfs) of the whole set and the pdfs of the samples in the different calibration sets. The msd values decreased as the sample set size increased. Despite this, the differences between the calibration sets in terms of their msd were marginal beyond 180 samples. Therefore, we used this number of samples as the optimal set size for the calibration of the vis-NIR models of the target soil properties. 

```r
## Optimal sample size for calibrating vis-NIR models
ocss <- 180
```

Now that we have the optimal set size, we select the very final calibration set using the conditioned Latin hypercube (cLH) algorithm. For this purpose, we sample 10 different sets with the inferred optimal size and selected the one with the minimum msd as the final calibration set.

```r
## Number of points in the density distribution
nxdens <- 2000

## create the matrix where the results will be stored
sc.dens <- data.frame(seq(min.sc[1], max.sc[1], length = nxdens), 
                      matrix(NA, nxdens, length(min.sc)))
colnames(sc.dens) <- c("x", paste("densc", 1:length(min.sc), sep = "")) 

## Estimate the density distribution of each component
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
for(i in 1:reps){
cat("------ Repetition", i, "------\n")
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

## Best group (obtained at iteration number 6)
calibration.idx <- calidx[[6]]
cal_smpls <- as.character(data$ID[calibration.idx])

## Then split the data of the candidate samples for calibration into 
## calibration (train) and prediction (pred) sets
train <- data[as.character(data$ID) %in%  cal_smpls, ]
pred <- data[!(as.character(data$ID) %in%  c(cal_smpls)),]
```

Now we have three basic subsets:
<ul>
<li> `train`: A calibration subset with an optimized size of 180 samples.</li>
<li> `pred`: A prediction subset of 503 samples. This subset will be used to predict the soil properties with the vis-NIR models that will be calibrated with `train`.</li>
<li> `valida`: A validation subset of 227 samples (corresponding to 115 sampling locations). This subset will be used to validate/assess the predictive performance of both, the vis-NIR models and the spatial models of soil proerties.</li></li>
</ul>

```r
dim(train)
dim(pred)
dim(valida)

train
pred
valida

## Save the results
save.image("sampling.RData")
``` 


## 2. Develop vis-NIR models to predict soil properties 
In this section, the idea is to use the `train` subset to develop vis-NIR models that allow us to predict the soil properties (sand, silt and clay contents as well as excangeable Ca) in the `pred` subset. Finally, for each soil property and for each layer, we will pool the reference wet-chemistry values (i.e. the values of the samples in `train`) and the vis-NIR predicted values (i.e. the values predicted for the samples in `pred`). Hereafter, we will refer to this resulting subset as vis-NIR augmented dataset. 
Here we will aos validate the predictive accuracy of the vis-NIR models using the `valida` subset.
<p>First, we will add a variable called layer that represents the depth at wich the sample was collected.</p>


```r
train$layer <- as.factor(substr(train$ID, start = 1, stop = 1))
pred$layer <- as.factor(substr(pred$ID, start = 1, stop = 1))
valida$layer <- as.factor(substr(valida$ID, start = 1, stop = 1))
``` 

Sand, silt and clay contents are reported as proportions, which sum up to 100%. However, predictive models built for each of these fractions do not guarantee that their individual predictions will sum up to 100%. To avoid this compositional constraint, the particle size data (clay,silt and sand) of both layers will be transformed and modeled based on the additive log-ratio (alr) transformation. To ilustrate how this transofrmation is applied we provide the following code example:

```r
## In this case we use the sand content as denominator to transform 
## the silt and clay content values
alr_clay <- log(train$Clay/train$Sand)
alr_silt <- log(train$Silt/train$Sand)

## For back-transformations we use:
clay <- 100 * exp(alr_clay)/(1 + exp(alr_clay) + exp(alr_silt))
silt <- 100 * exp(alr_silt)/(1 + exp(alr_clay) + exp(alr_silt))
sand <- 100 * 1/(1 + exp(alr_clay) + exp(alr_silt))

## As this piece of code is just to ilustrate how the alr transformation is applied
## we remove the objects that we just created
rm(list = c("alr_clay", "alr_silt", "clay", "silt", "sand"))

``` 
Each vis-NIR model will be developed using a memory-based learning (MBL) algorithm. This is, for each new sample (requiring prediction of a given soil property), a local model is calibrated using only its most similar samples (nearest neighbours) found in the calibration set. The most similar samples were selected based on similarity/dissimilarity among the spectra. Although MBL is commonly used to model large and complex vis-NIR datasets, it can also be used for modelling small vis-NIR datasets. In this case, rather than overcoming the typical complexity of large vis-NIR, the MBL algorithm is used to optimize the prediction of each sample by removing calibration samples that lay far apart from the prediction sample in the vis-NIR space. These distant samples can be considered as outliers, which may harm the accuracy of the prediction. 
<p> To build our MBL algorith we use the `resemble` package. First we define some basic apscts of our algorithm using the `mblControl` fucntion:

```r
## Ditance metric to select the neighbors or to remove local outliers:
## 'mahalanobis distance computed in the 'pls' factors (this is given by the sm argument)

## The method to select the number of components or PLS factors for distance computation will
## be the optimized pc selection ("opc") given in the pcSelection argument. The maximum number of 
## components to test will be limited to 40 (also given in the pcSelection argument)

## The (internal) validation method will be the Nearest neighbor validation
ctrl <- mblControl(sm = "pls", 
                   pcSelection = list("opc", 40), 
                   valMethod = c("NNv"), 
                   returnDiss = TRUE,
                   scaled = FALSE, 
                   center = TRUE)

## For additional details please consult the documentation of the resemble package
```

Now we define a vector with the threshold distances to test. Note: the number of threshold
distances tested here is large, it could be reduced
```r
diss2test <- seq(0.3, 1.5, by = 0.05)
```

Now we define the minimum (120) and maximum number of neighbors to be included in the
model.
```r
# Example: If with a given threshold distance only 70 neighbors are
# selected, then the function will be fored to take the especified minimum
# number of neighbors to be included in the model. In this case 120. In
# this case we also set that the maximum number of neighbors that can be
# included can be all the samples in the calibration set.
kminmax <- c(120, nrow(train$spc))
```

Now we will optimize the threshold distance for neighbor selection in MBL modles for each soil property. We will also validate the MBL models. 
<p>Exchangeable Ca:</p>

```r
## Here the calibration subset (i.e. train) is provided through the 
## arguments Xr (spectra) and Yr (reference value of the soil property) 
## while the spectral data of the subset for which the predictions of the 
## soil property are required (i.e. valida) is provided through the argument Xu. 

## The local regression method used here is the weighted average pls which 
## is given in the argument method. This is a case of ensemble modelling 
## where the final predictions are a weighted average of the predictions 
## generated by multiple (and consecutive) PLS components. Here we use the 
## average of the predictions retrieved by 16 PLS models 
## (ranging from 5 to 20 PLS components). This is provided in the argument pls.c.

## Note that we are not using the data of the response variable of the subset valida to 
## find the optimal distance threshold. 
sbl_ca <- mbl(Yr = train$Ca, 
              Xr = train$spc,
              Xu = valida$spc,
              mblCtrl = ctrl,
              group = substr(train$ID, 2, 100),
              dissUsage = "none",
              k.diss = diss2test, 
              k.range = kminmax,
              pls.c = c(5,20),
              method = "wapls1")
              
## Get the optimal threshold distance for the Ca predictions
idx.best.ca <- which.min(sbl_ca$nnValStats$st.rmse)
best.kdiss.ca <- sbl_ca$nnValStats$k.diss[idx.best.ca]

## Get the predicted values using the optimal threshold distance
pred_valCa <- getPredictions(sbl_ca)[,idx.best.ca]

## Now we can also validate the MBL predictions we just made
valR2Ca <- cor(valida$Ca, pred_valCa)^2 
valrmseCa <- (mean((valida$Ca - pred_valCa)^2))^0.5
valmeCa <- mean(valida$Ca - pred_valCa)

## Get the prediction variances (including the ones for each depth)
var.res.Ca_a <- var(valida[valida$layer=='A','Ca']
                    - pred_valCa[as.numeric(rownames(valida[valida$layer=='A',]))])
var.res.Ca_b <- var(valida[valida$layer=='B','Ca']
                    - pred_valCa[as.numeric(rownames(valida[valida$layer=='B',]))])
var.res.Ca <- var(valida$Ca - pred_valCa)

valR2Ca;valrmseCa;valmeCa

## Number of samples that are predicted by local regression (instead of
## using all the neighbors, i.e. global)
sum(sbl_ca$results[[idx.best.ca]]$k.org != nrow(train$spc))

## For additional details on the mbl function please check the documentation of the 
## resemble package
```

<p>Additive log-ratio transformation of silt:</p>

```r
## Run the MBL
sbl_alrsilt <- mbl(Yr = log(train$Silt/train$Sand), 
                   Xr = train$spc, 
                   Xu = valida$spc,
                   mblCtrl = ctrl,
                   group = substr(train$ID, 2, 100),
                   dissUsage = "none",
                   k.diss = diss2test, 
                   k.range = kminmax,
                   pls.c = c(5,20),
                   method = "wapls1")

## Get the optimal threshold distance for the log(Silt/Sand) predictions
idx.best.alrsilt <- which.min(sbl_alrsilt$nnValStats$st.rmse)
best.kdiss.alrsilt <- sbl_alrsilt$nnValStats$k.diss[idx.best.alrsilt]

## Get the predicted values using the optimal threshold distance
pred_val_alrsilt <- getPredictions(sbl_alrsilt)[,idx.best.alrsilt]
rownames(valida) <- seq(1,nrow(valida), by=1)

## Get the prediction variances (including the ones for each depth)
var.res.alr_silt_a <- var(valida[valida$layer=='A','alr_silt']
                          - pred_val_alrsilt[as.numeric(rownames(valida[valida$layer=='A',]))])
var.res.alr_silt_b <- var(valida[valida$layer=='B','alr_silt']
                          - pred_val_alrsilt[as.numeric(rownames(valida[valida$layer=='B',]))])
var.res.alrsilt <- var(valida$alr_silt - pred_val_alrsilt)

## Number of samples that are predicted by local regression (instead of 
## using all the neighbors, i.e. global)
sum(sbl_alrsilt$results[[idx.best.alrsilt ]]$k.org != nrow(train$spc))
```
<p>Additive log-ratio transformation of clay:</p>

```r
## Run the MBL
sbl_alrclay <- mbl(Yr = log(train$Clay/train$Sand), 
                   Xr = train$spc, 
                   Xu = valida$spc,
                   mblCtrl = ctrl,
                   group = substr(train$ID, 2, 100),
                   dissUsage = "none",
                   k.diss = diss2test, 
                   k.range = kminmax,
                   pls.c = c(5,20),
                   method = "wapls1")

## Get the optimal threshold distance for the log(Clay/Sand) predictions
idx.best.alrclay <- which.min(sbl_alrclay$nnValStats$st.rmse)
best.kdiss.alrclay <- sbl_alrclay$nnValStats$k.diss[idx.best.alrclay]

## Get the predicted values using the optimal threshold distance
pred_val_alrclay <- getPredictions(sbl_alrclay)[,idx.best.alrclay]

## Get the prediction variances (including the ones for each depth)
var.res.alr_clay_a <- var(valida[valida$layer=='A','alr_clay'] 
                          - pred_val_alrclay[as.numeric(rownames(valida[valida$layer=='A',]))])
var.res.alr_clay_b <- var(valida[valida$layer=='B','alr_clay']
                          - pred_val_alrclay[as.numeric(rownames(valida[valida$layer=='B',]))])
var.res.alrclay <- var(valida$alr_clay - pred_val_alrclay)

## Number of samples that are predicted by local regression (instead of 
## using all the neighbors, i.e. global)
sum(sbl_alrclay$results[[idx.best.alrclay]]$k.org != nrow(train$spc))
```

Back-transfrom to sand, silt and clay contents the predictions omade for the additive log-ratio transformations of the particle size. 

```r
valSilt_alr <- 100 * exp(pred_val_alrsilt)/(1+exp(pred_val_alrclay)+exp(pred_val_alrsilt))
valClay_alr <- 100 * exp(pred_val_alrclay)/(1+exp(pred_val_alrclay)+exp(pred_val_alrsilt))
valSand_alr <- 100 * 1/(1+exp(pred_val_alrclay)+exp(pred_val_alrsilt))
```

Compute the validation R^2^ and root mean squared error (RMSE) of the predicted sand, silt and clay contents

```r
valR2Silt_alr <- (cor(valida$Silt, valSilt_alr))^2
valR2Sand_alr <- (cor(valida$Sand, valSand_alr))^2
valR2Clay_alr <- (cor(valida$Clay, valClay_alr))^2

valrmseSilt_alr <- mean((valida$Silt - valSilt_alr)^2)^0.5
valrmseSand_alr <- mean((valida$Sand - valSand_alr)^2)^0.5
valrmseClay_alr <- mean((valida$Clay - valClay_alr)^2)^0.5

valmeSilt_alr <- mean(valida$Silt - valSilt_alr)
valmeSand_alr <- mean(valida$Sand - valSand_alr)
valmeClay_alr <- mean(valida$Clay - valClay_alr)
``` 
Now, we will use the optimal threshold distances to perform the MBL predictions of each soil property on the `pred` subset

```r
## Ca in the prediction subset
sbl_ca_pred <- mbl(Yr = train$Ca, 
                   Xr = train$spc, 
                   Xu = pred$spc,
                   mblCtrl = ctrl,
                   group = substr(train$ID, 2, 100),
                   dissUsage = "none",
                   k.diss = best.kdiss.ca, 
                   k.range = c(120, 180),
                   pls.c = c(5,20),
                   method = "wapls1")

## Get the predicted values
pred_predCa <- getPredictions(sbl_ca_pred)[,1]

## Number of samples that are predicted by local regression (instead of 
## using all the neighbors, i.e. global)
sum(sbl_ca_pred$results[[1]]$k.org != nrow(train$spc))

# Predict alr Silt in the prediction set based on the optimized threshold
# distance for neighbor selection
sbl_alrsilt_pred <- mbl(Yr = log(train$Silt/train$Sand), 
                        Xr = train$spc, 
                        Xu = pred$spc,
                        mblCtrl = ctrl,
                        group = substr(train$ID, 2, 100),
                        dissUsage = "none",
                        k.diss = best.kdiss.alrsilt, 
                        k.range = c(120, 180),
                        pls.c = c(5,20),
                        method = "wapls1")

## Get the predicted values
pred_pred_alrsilt <- getPredictions(sbl_alrsilt_pred)[,1]

## Number of samples that are predicted by local regression (instead 
# of using all the neighbors, i.e. global)
sum(sbl_alrsilt_pred$results[[1]]$k.org != nrow(train$spc))

# Predict alr Clay in the prediction set based on the optimized threshold
# distance for neighbor selection
sbl_alrclay_pred <- mbl(Yr = log(train$Clay/train$Sand), 
                        Xr = train$spc, 
                        Xu = pred$spc,
                        mblCtrl = ctrl,
                        group = substr(train$ID, 2, 100),
                        dissUsage = "none",
                        k.diss = best.kdiss.alrclay, 
                        k.range = c(120, 180),
                        pls.c = c(5,20),
                        method = "wapls1")

## Get the predicted values
pred_pred_alrclay <- getPredictions(sbl_alrclay_pred)[,1]

## Number of samples that are predicted by local regression (instead 
## of using all the neighbors, i.e. global)
sum(sbl_alrclay_pred$results[[1]]$k.org != nrow(train$spc))

## back transfrom the additive log-ratio transformations of the particle size. 
predSilt_alr <- 100 * exp(pred_pred_alrsilt)/(1+exp(pred_pred_alrclay)+
                                                exp(pred_pred_alrsilt))
predClay_alr <- 100 * exp(pred_pred_alrclay)/(1+exp(pred_pred_alrclay)+
                                                exp(pred_pred_alrsilt))
predSand_alr <- 100 * 1/(1+exp(pred_pred_alrclay)+exp(pred_pred_alrsilt))
```
Now let's store the results


```r
results_spec_modeling <- data.frame(Property = c("Sand", "Silt", "Clay", "Ca"), 
                                    R2 = rep(NA, 4), 
                                    RMSE = rep(NA, 4),
                                    ME = rep(NA, 4))

results_spec_modeling$R2 <- c(valR2Sand_alr, 
                              valR2Silt_alr, 
                              valR2Clay_alr, 
                              valR2Ca)
results_spec_modeling$RMSE <- c(valrmseSand_alr, 
                                valrmseSilt_alr, 
                                valrmseClay_alr, 
                                valrmseCa)
results_spec_modeling$ME <- c(valmeSand_alr, 
                              valmeSilt_alr, 
                              valmeClay_alr, 
                              valmeCa)

## Write a table of the results
write.table(results_spec_modeling, 
            file = "results_spec_modeling.txt", 
            row.names = FALSE, sep = "\t")
``` 
Plot the results
```r 
locsamples <- rbind(data.frame(ns = as.numeric(as.character(sbl_alrsilt_pred$results[[1]]$k)), 
                               Property = "alrsilt"),
                    data.frame(ns = as.numeric(as.character(sbl_alrclay_pred$results[[1]]$k)), 
                               Property = "alrclay"),
                    data.frame(ns = as.numeric(as.character(sbl_ca_pred$results[[1]]$k)), 
                               Property = "Ca"))


sum(sbl_ca_pred$results[[1]]$k.org == nrow(train$spc))/nrow(pred)
sum(sbl_alrsilt_pred$results[[1]]$k.org == nrow(train$spc))/nrow(pred)
sum(sbl_alrclay_pred$results[[1]]$k.org == nrow(train$spc))/nrow(pred)

ggplot(locsamples, aes(x = ns, fill = Property)) +
  geom_histogram(aes(y = ..density..), 
                 alpha = 0.5,
                 position = 'identity', 
                 binwidth=10)
```
Prepare the vis-NIR-augemnted databas for spatial modelling
```r
specaidaddb <- rbind(data.frame(train[,c("ID", 
                                         "POINT_X", 
                                         "POINT_Y", 
                                         "layer", 
                                         "Clay", 
                                         "Sand", 
                                         "Silt", 
                                         "Ca")], 
                                set = "train", 
                                specClay_alr = train$Clay, 
                                specSand_alr = train$Sand, 
                                specSilt_alr = train$Silt,
                                specCa = train$Ca), 
                     
                     data.frame(valida[,c("ID", 
                                          "POINT_X", 
                                          "POINT_Y", 
                                          "layer", 
                                          "Clay", 
                                          "Sand", 
                                          "Silt", 
                                          "Ca")], 
                                set = "validation", 
                                specClay_alr = valida$Clay, 
                                specSand_alr = valida$Sand, 
                                specSilt_alr = valida$Silt,
                                specCa = valida$Ca),  
                     
                     data.frame(pred[,c("ID", 
                                        "POINT_X", 
                                        "POINT_Y", 
                                        "layer", 
                                        "Clay", 
                                        "Sand", 
                                        "Silt", 
                                        "Ca")], 
                                set = "prediction", 
                                specClay_alr = predClay_alr, 
                                specSand_alr = predSand_alr, 
                                specSilt_alr = predSilt_alr,
                                specCa = pred_predCa))
                                
## Write table of the vis-NIR database
write.table(specaidaddb, file = "spec_aideddb.txt", row.names = F, sep = "\t")
```
Plot and save the additional data
```r
wav <- as.numeric(colnames(train$spc))
spca <- 1/exp(data$spc[substr(data$ID, 1, 1) == "A", ])
spcb <- 1/exp(data$spc[substr(data$ID, 1, 1) == "B", ])

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

write.table(round(cormat, 2), 
            file = "correlation_matrix.txt", 
            row.names = FALSE, 
            sep = "\t")

## Save the variance of the spectroscopic model error
save(var.res.Ca, 
     var.res.Ca_a, 
     var.res.Ca_b, 
     var.res.alr_clay_a,
     var.res.alr_clay_b,
     var.res.alrclay,
     var.res.alr_silt_a, 
     var.res.alr_silt_b, 
     var.res.alrsilt, 
     file = "var_spec_error.Rdata")
```

## 3. Spatial modeling of soil properties using the vis-NIR augemted data
```r
data.spec <- read.table("spec_aideddb.txt", 
                        header = TRUE, 
                        check.names = FALSE, 
                        sep ="\t")
```


