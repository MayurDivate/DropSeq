# dropseq data 
# genes and trascript violin plot 
# for batch processing 
# author@MayurDivate


# set path to dir where Human and Mouse dge.summary files 
# are located for all samples 

setwd("/path/summary")

hs <- list.files(pattern = ".*Human.*txt")
ms <- list.files(pattern = ".*Mouse.*txt")

files <- data.frame(hs,ms,stringsAsFactors = FALSE)

for( i in 1:nrow(files)){
  
  HumanFile <- files[i,1]
  MouseFile <- files[i,2]
  
  plotTitle <- gsub("_Human.dge.summary.txt","",HumanFile)
  plotGimg <- paste(plotTitle,"_genes.jpg",sep = "")
  plotTimg <- paste(plotTitle,"_transcripts.jpg",sep = "")
  
  dfH <- read.table(HumanFile, sep="\t", header=TRUE, stringsAsFactors = F)
  dfM <- read.table(MouseFile, sep="\t", header=TRUE, stringsAsFactors = F)
  
  colnames(dfH) <- c("CELL_BARCODE","HGenes","HTrancripts")
  colnames(dfM) <- c("CELL_BARCODE","MGenes","MTrancripts")
  
  dfHM <- merge(dfH,dfM,by="CELL_BARCODE")
  dfHM$totalG <- rowSums(dfHM[,c(2,4)])
  dfHM$totalT <- rowSums(dfHM[,c(3,5)])
  dfHM$RatioH <- dfHM$HTrancripts / dfHM$totalT

  dfHM$Species <- rep("Mixed",nrow(dfHM)) 
  dfHM$Species[dfHM$RatioH > 0.8] = "Human"
  dfHM$Species[dfHM$RatioH < 0.2 ] = "Mouse"
  
  dfHM$transcripts <- rep(0,nrow(dfHM))
  dfHM$genes <- rep(0,nrow(dfHM))
  
  dfHM[dfHM$Species == "Human",10] = dfHM[dfHM$Species == "Human",3]
  ht <- round(mean(dfHM[dfHM$Species == "Human",3]),2)
  dfHM[dfHM$Species == "Mouse",10] = dfHM[dfHM$Species == "Mouse",5]
  mt <- round(mean(dfHM[dfHM$Species == "Mouse",5]),2)
  
  dfHM[dfHM$Species == "Human",11] = dfHM[dfHM$Species == "Human",2]
  hg <- round(mean(dfHM[dfHM$Species == "Human",2]),2)
  dfHM[dfHM$Species == "Mouse",11] = dfHM[dfHM$Species == "Mouse",4]
  mg <- round(mean(dfHM[dfHM$Species == "Mouse",4]),2)
  dfV <- dfHM[dfHM$Species != "Mixed",c(9,10,11)] 
  
  c <- data.frame(sample =c(plotTitle) , ht ,mt,hg, mg)
  
  jpeg(plotGimg,width = 1000, height = 1000, quality = 100,res=200)

  plot <- ggplot(dfV,aes( Species, genes,colour = Species))
  plot <- plot + geom_violin()
  plot <- plot + labs(x = "Species", y= "Genes")
  plot <- plot + labs(title = plotTitle, subtitle = paste(plotTitle," Human = ",hg,", Mouse = ",mg,sep = ""))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5) )
  plot <- plot + theme(legend.position = c(.95, .95),
                       legend.justification = c("right", "top"),
                       legend.box.just = "right",
                       legend.margin = margin(6, 6, 6, 6))
  
  plot <- plot + geom_boxplot(width = 0.1,outlier.color = "black", fill="grey")
  
  options(scipen=10000)
  
  print(plot)
  
  dev.off()
  
  
  
  jpeg(plotTimg,width = 1000, height = 1000, quality = 100,res=200)
  
  plot <- ggplot(dfV,aes( Species, transcripts,colour = Species))
  plot <- plot + geom_violin()
  plot <- plot + labs(x = "Species", y= "Transcripts")
  plot <- plot + labs(title = plotTitle, subtitle = paste(plotTitle," Human = ",ht,", Mouse = ",mt,sep = ""))
  plot <- plot + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  plot <- plot + theme(legend.position = c(.95, .95),
                       legend.justification = c("right", "top"),
                       legend.box.just = "right",
                       legend.margin = margin(6, 6, 6, 6))
  
  plot <- plot + geom_boxplot(width = 0.1,outlier.color = "black", fill="grey")
  
  options(scipen=10000)
  
  print(plot)
  
  dev.off()

}


