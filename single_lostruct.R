setwd("~/3.SV/13.localPCA")
args = commandArgs(trailingOnly=TRUE)
#setwd("/Users/rainlam")
library(ggplot2)
library(tidyverse)
library(dplyr)


###########################################################
################## 1. local PCA ###########################
###########################################################
#pop <- "Aripo" #
pop <- args[1]
#chr <- "LG4"   #
chr <- args[2]
window_size <- 100
min_windows <- 3
n_permutations <- 1000
k_kept <- 40

#lostruct
devtools::load_all("/home/linyuying/bin/local_pca/lostruct/R")
library(lostruct)

bcf.file <- paste(pop, "/", pop, ".", chr, ".vcf.gz", sep="")
# read VCF file
sites <- vcf_positions(bcf.file)
#separate vcf into windows
win.fn.snp <- vcf_windower(bcf.file, size=window_size, type="snp", sites=sites)
#PCA in windows
system.time(snp.pca <- eigen_windows(win.fn.snp, k=2))
# PC distance
system.time(pcdist <- pc_dist(snp.pca, mc.cores=10))
# remove missing data
pcdist_na <- which(is.na(pcdist), TRUE)
na.inds <- is.na(pcdist[,1]) 
if (sum(na.inds) == length(na.inds)){
  na.inds <- is.na(pcdist[,2]) 
}
mds <- cmdscale( pcdist[!na.inds,!na.inds], eig=TRUE, k=k_kept)

# mds for each PC
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.regions <- region(win.fn.snp)()
win.mids <- (win.regions$start+win.regions$end)/2
win.regions$n <- 1:nrow(win.regions)
win.regions <- win.regions[!na.inds,]
win.regions %>% mutate(mid = (start + end) / 2) ->  win.regions
# Add the columns for all the MDS coordinates
for (k in 1:k_kept){
  str_pad(k, 2, pad = "0")
  
  name = paste("mds",str_pad(k, 2, pad = "0"),sep="")
  win.regions$tmp <- "NA"
  win.regions <- win.regions %>% rename(!!name := tmp) 
}

#Add the MDS coordinates to each window.
for (i in 1:k_kept){
  j = i + 5
  win.regions[,j] <- mds.coords[,i]
}

#save RDS
saveRDS(win.regions, file = paste(pop, ".", chr, ".win", window_size, ".lostruct.windows.rds", sep="")) 
saveRDS(mds, file = paste(pop, ".", chr, ".win",  window_size, '.mds', sep=""))


###########################################################
#################### corner ###############################
###########################################################
mincirc <- enclosing_circle( mds.coords[,1:2] )
mds.corners <- corners( mds.coords[,1:2], prop=.05 )
corner.cols <- c("red","blue","purple")
ccols <- rep("black",nrow(mds.coords))
for (k in 1:ncol(mds.corners)) {
  ccols[ mds.corners[,k] ] <- corner.cols[k]
}
pdf(paste(pop, ".", chr, ".win", window_size,".corner.pdf", sep=""))
plot( mds.coords[,1:2], pch=20, col=adjustcolor(ccols,0.75), 
      xlab="MDS coordinate 1", ylab="MDS coordinate 2", 
      xlim=mincirc$ctr[1]+c(-1,1)*mincirc$rad,
      ylim=mincirc$ctr[2]+c(-1,1)*mincirc$rad )
plot_circle( mincirc, col='red' )
points( mincirc$three, col='red', cex=2 )
points( mds.coords[mincirc$index,], col='red', cex=1.5 )
dev.off()
# new win.regions
win.regions$cclo <- ccols

pdf(paste(pop, ".", chr, ".win", window_size,".lostruct.MDSplots.pdf",sep=""), height=16,width=25)
print(
  win.regions %>%
    gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
    ggplot(.,aes(x=mid,y=value, colour=cclo)) + geom_point() +
    scale_color_manual(values = c("black" = "black", "purple" = "purple","blue"="blue","red"="red")) + 
    facet_grid(mds~.,scales = "free") +
    theme_bw()
)

dev.off()


##########################################################
##################### 3. plotting ########################
##########################################################
for (k in 6:15) {
  pdf(paste(pop, chr, "win", window_size,"MDS", k - 5, "pdf", sep ="." ))
  plot( win.regions$mid/1e6, win.regions[,k], pch=20, 
        xlab="Position (Mb)", ylab=paste("MDS coordinate", k - 5),
        col=adjustcolor(ccols,0.75) )
  dev.off()
}

