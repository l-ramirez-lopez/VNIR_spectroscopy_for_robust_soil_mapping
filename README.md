# R code for 'Robust soil mapping at the farm scale with vis-NIR spectroscopy'

_Leo Ramirez-Lopez and Alex Wadoux_

_::::05.06.2019::::_


```{r setup, include=FALSE}
library(leaflet)
library(dplyr)
```

Use the leaflet map below to explore the actual Maunga Whau volcano in Auckland, NZ. 

```{r}
leaflet() %>%
  setView(lng=174.764, lat=-36.877, zoom = 16) %>% 
  addTiles() %>%
  addMarkers(lng=174.764, lat=-36.877, popup="Maunga Whau") 
```


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## GitHub Documents

This is an R Markdown format used for publishing markdown documents to GitHub. When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated.

## Including Code

You can include R code in the document as follows:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=TRUE}
plot(pressure)
```



In the spirit of reproducible research, here we share the data and the computational code
used for carrying out the analyses presented in our paper ["Robust soil mapping at the farm scale with vis-NIR spectroscopy"](https://onlinelibrary.wiley.com/doi/10.1111/ejss.12752) which is by the way open access. 

### paper summary
Sustainable agriculture practices are often hampered by the prohibitive costs associated with the generation of fine-resolution soil maps. Recently, several papers have been published highlighting how visible and near infrared(vis-NIR) reflectance spectroscopy may offer an alternative to address this problem by increasing the density of soil sampling and by reducing the number of conventional laboratory analyses needed. However, for farm-scale soil mapping, previous studies rarely focused on sample optimization for the calibration of vis-NIR models or on robust modelling of the spatial variation of soil properties predicted by vis-NIR spectroscopy. In the present study, we used soil vis-NIR spectroscopy models optimized in terms of both number of calibration samples and accuracy for high-resolution robust farm-scale soil mapping and addressed some of the most common pitfalls identified in previous research. We collected 910 samples from 458 locations at two depths (A, 0-0.20 m; B, 0.80-1.0 m) in the state of S�o Paulo, Brazil. All soil samples were analysed by conventional methods and scanned in the vis-NIR spectral range. With the vis-NIR spectra only, we inferred statistically the optimal set size and the best samples with which to calibrate vis-NIR models. The calibrated vis-NIR models were validated and used to predict soil properties for the rest of the samples. The prediction error of the spectroscopic model was propagated through the spatial analysis, in which robust block kriging was used to predict particle-size fractions and exchangeable calcium content for each depth. The results indicated that statistical selection of the calibration samples based on vis-NIR spectra considerably decreased the need for conventional chemical analysis for a given level of mapping accuracy.

### In case of comments/questions/issues 
In case you have any question comments/questions/issue related to the code presented here please go to:  
[https://github.com/l-ramirez-lopez/VNIR_spectroscopy_for_robust_soil_mapping/issues](https://github.com/l-ramirez-lopez/VNIR_spectroscopy_for_robust_soil_mapping/issues) and create a new issue.

### For citation or details please refer to: 

Ramirez-Lopez, L., Wadoux, A. C., Franceschini, M. H. D., Terra, F. S., Marques, K. P. P., Say�o, V. M., & Dematt�, J. A. M. (2019). [Robust soil mapping at the farm scale with vis-NIR spectroscopy](https://onlinelibrary.wiley.com/doi/10.1111/ejss.12752). European Journal of Soil Science.

### Notes
* _Reproducibility_: slight discrepancies between the results of reoported in the paper and the results you get in your `R` console might be expected in case you use different `R` package versions and different random number generators.  

* _To advanced R users_: In order to make the `R` code more readable/interpretable we have opted for using for using for loops instead of vectorized functions. There are several aspects of the code that can be improved to reach better computational efficiency.



### code organization
* 1. Sampling: here we cover the section covers the section of the paper 'Selection of calibration and validation sets'
* 2. 






# 1. Sampling


In the code below is used to select the calibration and validation sets

------------------------------------------------------------------------

Install the required packages

``` r
requiredpackages <- c("resemble",    # version 1.2.2
                      "prospectr",   # version 0.1.3
                      "clhs",        # version 0.5-6
                      "matrixStats", # version 0.51.0
                      "doParallel",  # version 1.0.10
                      "ggplot2"      # version 2.2.1
)

    
toinstall <- requiredpackages[!requiredpackages %in% rownames(installed.packages())]

if(length(toinstall) > 0){
  install.packages(toinstall)
}

lapply(requiredpackages, 
       FUN = library, 
       character.only = TRUE)
```

Define the working directory

``` r
workingd <- "/myworkingdirectory"
setwd(workingd)
```

Organize the data
=================

Read the data set

``` r
# In case the data is already saved locally
data <- readRDS("SoilNIRSaoPaulo.rds")
```
or

``` r
# In case you want to read the data directly from the repository
nirfile <- file("https://github.com/l-ramirez-lopez/VNIR_spectroscopy_for_robust_soil_mapping/raw/master/SoilNIRSaoPaulo.rds")

data <- readRDS(nirfile)
```

Apply standard normal variate

``` r
data$spc_snv <- standardNormalVariate(data$spc)
```

Remove the validation samples from the analysis

``` r
data <- data[data$set == "cal_candidate",]
```

Perform a principal component analysis
======================================

Compress the data

``` r
pcaall <- orthoProjection(Xr = data$spc_snv, 
                          X2 = NULL, 
                          Yr = NULL, 
                          method = "pca", 
                          pcSelection = list("cumvar", 0.99), 
                          center = TRUE, 
                          scaled = FALSE)
```

Standardize the scores

``` r
pcaall$scores.std <- sweep(pcaall$scores, MARGIN = 2, STATS = pcaall$sc.sdv, FUN = "/")
```

Compute the maximum and minimum of each score for the limits of the density estimations

``` r
max.sc <- colMins(pcaall$scores.std)
min.sc <- colMaxs(pcaall$scores.std)
```

Compute the mean and standard deviation of each score for the comparisons with the samples

``` r
mean.sc <- colMeans(pcaall$scores.std)
sd.sc <- colSds(pcaall$scores.std)
```

Number of points in the density distribution

``` r
nxdens <- 500
```

Matrix where the density values will be stored

``` r
ix <- 1
sc.dens <- data.frame(seq(min.sc[ix], max.sc[ix], length = nxdens), 
                      matrix(NA, nxdens, length(min.sc)))
colnames(sc.dens) <- c("x", paste("densc", 1:length(min.sc), sep = ""))
```

Estimate the density distribution of each component

``` r
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

Create different sample set size

``` r
css <- seq(10, 400, by = 10)
```

Sample with conditioned latin hypercube sampling (clhs)

``` r
# Sampling repetitions
sreps <- 10

results.clhs <-  data.frame(css = css, 
                            msd = rep(NA, length(css)),
                            mndiff = rep(NA, length(css)),
                            sddiff = rep(NA, length(css)))

# these two nested loops will take a while (2 nested loops:
# one for testing 40 different calibration set sizes and the other
# to repeat 10 times the whole procedure)

# in this case we have use loops because of 
# code interpretability reasons
# R tip:
# this part of the code could be re-written to speed up computations
# e.g. the operations here can be vectorized by using the family of
# apply functions and also by using parallel computing


for(k in 1:sreps){
  results.clhs[,-1] <- NA
  fn <- paste("6pcs_resultsclhs_rep", k,".txt", sep = "")
  
  if(fn %in% list.files())
  {
    results.clhs <- read.table(fn, header = T, sep = "\t")
  }
  
  iter.p <- 1 + sum(rowSums(!is.na(results.clhs)) == ncol(results.clhs))
  
  for(i in iter.p:length(css)){
    set.seed(k)
    i.calidx <- clhs(as.data.frame(pcaall$scores.std),
                     size = css[i], iter = 10000)
    
    j.sc.dens <- msd.sc  <- sc.dens 
    for(j in 1:length(min.sc)){
      
      # use the same bandwidth (bw) as in the whole set of candidates
      j.sc.dens[,j+1] <- density(pcaall$scores.std[i.calidx,j], 
                                 bw = d.bandwidths[j], 
                                 n = nxdens, from = min.sc[j], to = max.sc[j], 
                                 kernel = "gaussian")$y
    }
    results.clhs$msd[i] <- mean(colMeans((j.sc.dens[,-1] - sc.dens[,-1])^2, na.rm = T))
    results.clhs$mndiff[i] <- mean(abs(colMeans(pcaall$scores.std[i.calidx,])))
    results.clhs$sddiff[i] <- mean(abs(colSds(pcaall$scores.std[i.calidx,]) - 1))
    
    write.table(results.clhs, 
                file = fn,
                row.names = FALSE, sep = "\t")
    
    print(results.clhs[1:i,])
  }
}

nmsrepsclhs <- paste("6pcs_resultsclhs_rep", 1:sreps, ".txt", sep = "")
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
```

Plot the results using ggplot

``` r
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
```

Save the plot

``` r
pdf(file = "calibration_set_size.pdf", width = 8, height = 6)
p.tmp + geom_ribbon(aes(ymin=sd_lower, ymax=sd_upper, x=css, colour = "bands"), alpha = 0.2) +
  scale_colour_manual(name = '', values = c("bands" = NA, "msd" = "black"))
dev.off()
```

Save the results

``` r
final.clhs <- read.table("6pcs_final.clhs.txt", header = TRUE, check.names = FALSE, sep ="\t")
save.image("sampling.RData")
```

