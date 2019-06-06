1. Sampling
================

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
workingd <- "myworkingdirectory/"
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

# these two nested loops will take a while
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
