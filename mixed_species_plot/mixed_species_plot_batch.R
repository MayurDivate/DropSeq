# for batch processing 
# author@MayurDivate
 

# set path to dir where Human and Mouse dge.summary files 
# are located for all samples 
  
setwd("/path/summary")

hs <- list.files(pattern = ".*Human.*txt")
ms <- list.files(pattern = ".*Mouse.*txt")

files <- data.frame(hs,ms,stringsAsFactors = FALSE)
head(files)

for( i in 1:nrow(files)){
  
HumanFile <- files[i,1]
MouseFile <- files[i,2]

print(HumanFile)
print(MouseFile)

plotTitle <- gsub("_Human.dge.summary.txt","",HumanFile)
plotimg <- paste(plotTitle,"_mixed_species.jpg",sep = "")
print(plotimg)

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

cellCounts <- as.data.frame(summary(factor(dfHM$Species)))

  if(length(cellCounts[rownames(cellCounts) == "Human",]) > 0){
    hg <- paste("Human (", cellCounts["Human",],")")
    print(hg)
  }

  if(length(cellCounts[rownames(cellCounts) == "Mouse",]) > 0){
    mm <- paste("Mouse (", cellCounts["Mouse",],")")
    print(mm)
  }

  if(length(cellCounts[rownames(cellCounts) == "Mixed",]) > 0){
    mx <- paste("Mixed (", cellCounts["Mixed",],")")
    print(mx)
  }

dfHM$Species[dfHM$Species == "Human"] = hg 
dfHM$Species[dfHM$Species == "Mouse"] = mm 
dfHM$Species[dfHM$Species == "Mixed"] = mx 


jpeg(plotimg,width = 1000, height = 1000, quality = 100,res=200)
plot <- ggplot(dfHM,aes( HTrancripts, MTrancripts, colour = Species))
plot <- plot + geom_point(size = 1)
plot <- plot + labs(x = "Human transcripts", y= "Mouse transcripts")
plot <- plot + labs(title = plotTitle)
plot <- plot + theme(plot.title = element_text(hjust = 0.5))
plot <- plot + theme(legend.position = c(.95, .95),
                     legend.justification = c("right", "top"),
                     legend.box.just = "right",
                     legend.margin = margin(6, 6, 6, 6))

options(scipen=10000)

print(plot)

dev.off()

}

# done 