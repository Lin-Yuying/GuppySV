setwd("~/localPCA/")
setwd("/Users/rainlam")
library(ggplot2)
library(tidyverse)
pops <-c("Aripo_Quare", "Aripo_Yarra", "Quare_Yarra")
chr <- paste('LG', seq(1,23), sep = '')
window.size <- 100

devtools::load_all("~/localPCA/local_pca/lostruct/R")

for (pop in pops) {
  genome.chr <- paste(pop, chr, sep='.') #"Quare_Yarra.LG1" 
  win.regions <- readRDS(paste(genome.chr, window.size, "lostruct.windows.rds", sep='.'))
  mds <- readRDS(paste(genome.chr, window.size,'mds', sep='.'))
  mds.coords <- mds$points
  colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
  win.mids <- (win.regions$start+win.regions$end)/2
  mincirc <- enclosing_circle( mds.coords[,1:2] )
  mds.corners <- corners( mds.coords[,1:2], prop=.05 )
  corner.cols <- c("red","blue","purple")
  ccols <- rep("black",nrow(mds.coords))
  for (k in 1:ncol(mds.corners)) {
    ccols[ mds.corners[,k] ] <- corner.cols[k]
  }
  # MD1 and MD2 plot
  pdf(paste(genome.chr, window.size, 'MDs.pdf', sep="."))
  plot( mds.coords[,1:2], pch=20, col=adjustcolor(ccols,0.75), 
        xlab="MDS coordinate 1", ylab="MDS coordinate 2", 
        xlim=mincirc$ctr[1]+c(-1,1)*mincirc$rad,
        ylim=mincirc$ctr[2]+c(-1,1)*mincirc$rad,
        main = 'Yarra.LG1')
  plot_circle( mincirc, col='red' )
  points( mincirc$three, col='red', cex=2 )
  points( mds.coords[mincirc$index,], col='red', cex=1.5 )
  dev.off()
  
  # MDs color plots 
  pdf(filename = paste(genome.chr, window.size,"MDs.color.pdf",sep="."), height = 2000, width = 1000, res = 100)
  print(
    win.regions$cclor = ccols
    win.regions %>%
      gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
      ggplot(.,aes(x=mid,y=value)) + geom_point(aes(color = cclor), show.legend = FALSE) + facet_grid(mds~.,scales = "free") +
      scale_color_manual(values = c("black" = "black", "blue" = "blue", "purple" = 'purple', "red" = "red")) +
      theme_bw()
  )
  dev.off()
  # MDs plot without color
  pdf(filename = paste(genome.chr, window.size, 'MDs.pdf', sep="."), height = 2000, width = 1000, res = 100)
  print(
    win.regions %>%
      gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
      ggplot(.,aes(x=mid,y=value)) + geom_point() + facet_grid(mds~.,scales = "free") +
      theme_bw()
  )
  dev.off()
  
}


win.regions <- readRDS("Aripo.LG4.100.lostruct.windows.rds")
mds <- readRDS('Aripo.LG4.100.mds')
mds.coords <- mds$points
colnames(mds.coords) <- paste("MDS coordinate", 1:ncol(mds.coords))
win.mids <- (win.regions$start+win.regions$end)/2
mincirc <- enclosing_circle( mds.coords[,1:2] )
mds.corners <- corners( mds.coords[,1:2], prop=.05 )
corner.cols <- c("red","blue","purple")
ccols <- rep("black",nrow(mds.coords))
for (k in 1:ncol(mds.corners)) {
  ccols[ mds.corners[,k] ] <- corner.cols[k]
}

pdf("Aripo.LG4.MDs.pdf")
plot( mds.coords[,1:2], pch=20, col=adjustcolor(ccols,0.75), 
      xlab="MDS coordinate 1", ylab="MDS coordinate 2", 
      xlim=mincirc$ctr[1]+c(-1,1)*mincirc$rad,
      ylim=mincirc$ctr[2]+c(-1,1)*mincirc$rad,
      main = 'Aripo.LG4')
plot_circle( mincirc, col='red' )
points( mincirc$three, col='red', cex=2 )
points( mds.coords[mincirc$index,], col='red', cex=1.5 )
dev.off()
# MDs plots 
pdf("Aripo.LG4.mds.color.pdf", height = 16, width = 25)
win.regions$cclor = ccols
print(
  win.regions %>%
    gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
    ggplot(.,aes(x=mid,y=value)) + geom_point(aes(color = cclor), show.legend = FALSE) + facet_grid(mds~.,scales = "free") +
    scale_color_manual(values = c("black" = "black", "blue" = "blue", "purple" = 'purple', "red" = "red")) +
    theme_bw()
)

dev.off()

pdf("Aripo.LG4.mds.pdf", height = 16, width = 25)
print(
  win.regions %>%
    gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
    ggplot(.,aes(x=mid,y=value)) + geom_point() + facet_grid(mds~.,scales = "free") +
    theme_bw()
)
dev.off()

# plot 10 MDs
win.regions$cclor = ccols
win.regions %>%
  gather(., mds, value, colnames(win.regions)[6:(ncol(win.regions)-30)]) %>% 
  ggplot(.,aes(x=mid,y=value)) + geom_point(aes(color = cclor), show.legend = FALSE) + facet_grid(mds~.,scales = "free") +
  scale_color_manual(values = c("black" = "black", "blue" = "blue", "purple" = 'purple', "red" = "red")) +
  theme_bw()
  
# single plot
ggplot(data = df, aes(x = win.mids/1e6, y = as.numeric(MDs1))) + 
  geom_point(aes(color = clr), show.legend = FALSE) +
  xlab("LG1 / MB") + ylab("mds1") +
  scale_color_manual(values = c("black" = "black", "blue" = "blue", "purple" = 'purple', "red" = "red")) +
  theme_bw()
  
corner.npc <- 4
corner.regions <- lapply( 1:ncol(mds.corners), function (k) {
  regions[ mds.corners[,k],]
} )
corner.covmats <- lapply( 1:ncol(mds.corners), function (k) {
  reg <- regions[ mds.corners[,k], ]
  qfun <- multi_vcf_query_fn( chrom.list=chroms, file=bcf.files, regions=reg )
  running_cov(qfun,1:nrow(reg))
} )
corner.pca <- lapply( corner.covmats, function (covmat) {
  cov_pca(covmat=covmat,k=corner.npc,w=opt$weights)
} )  
