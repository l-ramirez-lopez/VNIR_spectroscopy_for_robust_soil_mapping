Annexe 1
================

========================================================

Refers to section Selecting the calibration and validation sets

------------------------------------------------------------------------

Install the required packages

``` r
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
options(warn=0)
```

Define the working directory

``` r
workingd <- "your/working/directory/"
setwd(workingd)
```

Organize the data
=================

Read the data

``` r
data <- read.table("SoilNIRSaoPaulo.txt", header = TRUE, check.names = FALSE, sep ="\t")
```

Minimum wavelength in the dataset

``` r
mwav <- 502
spc <- log(1/data[,which(colnames(data) == mwav):ncol(data)])
data <- data[,-c(which(colnames(data) == mwav):ncol(data))]
data$spc <- spc
```

Apply standard normal variate

``` r
data$spc_snv <- standardNormalVariate(data$spc)
rm(spc)
```

Remove the validation samples from the analysis

``` r
data <- data[data$set == "cal_candidate",]
```

Perform a principal component analysis
======================================

Compress the data

``` r
pcaall <- orthoProjection(Xr = data$spc_snv, X2 = NULL, 
                          Yr = NULL, 
                          method = "pca", pcSelection = list("cumvar", 0.99), 
                          center = TRUE, scaled = FALSE, cores = 1)
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
results.clhs <-  data.frame(css = css, 
                            msd = rep(NA, length(css)),
                            mndiff = rep(NA, length(css)),
                            sddiff = rep(NA, length(css)))

for(k in 1:10){
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

nmsrepsclhs <- paste("6pcs_resultsclhs_rep", 1:10, ".txt", sep = "")
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
