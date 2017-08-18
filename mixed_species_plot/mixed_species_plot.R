# Before using this script you should finsh analsysis using dropseq tools
# dropseq tools : http://mccarrolllab.com/dropseq/
# Use dropseq cookbook to analyze data using dropseq tools
# perform digital gene expression analysis 
# get 2 dge files one for Human and one for Mouse
# author@MayurDivate

# human digital gene expression file
HumanFile <- "SRR1873277_Human.dge.summary.txt"
# mouse digital gene expression file
MouseFile <- "SRR1873277_Mouse.dge.summary.txt"


plotTitle <- gsub("_Human.dge.summary.txt","",HumanFile)
plotimg <- paste(plotTitle,"_mixed_species.jpg",sep = "")

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

# plotting data 

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

# done