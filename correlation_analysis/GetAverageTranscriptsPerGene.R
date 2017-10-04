# this script will give the matrix of Avergae genes per million
# for all the samples in specified folder 
# author@MayurDivate
# 


# specifiy the folder path
setwd("~/Work/0_GitHub/DropSeq/correlation_analysis/sample_data/")

# output file name 
outputFile <- "AverageTrancriptsPerGene.txt"

files <- list.files(pattern = "*.gz")
sampleDGEFile <- files[1]
sampleDGEFile
sample <- "TEST"

sampleDF <- read.table(sampleDGEFile,header = T,sep = "\t",stringsAsFactors = F)
rownames(sampleDF) <- sampleDF[,1]
sampleDF <- sampleDF[,c(2:ncol(sampleDF))]
sampleDF$mean  <- rowMeans(sampleDF)
sampleDF$genes <- rownames(sampleDF)
gCol <- ncol(sampleDF)
mCol <- gCol - 1

head(sampleDF)

resultDF <- sampleDF[,c(gCol,mCol)]
rownames(resultDF) <- c(1:nrow(resultDF))
colnames(resultDF) <- c("genes",sample)

for(sampleDGEFile in files){
  sample <- gsub(".dge.txt.gz","",sampleDGEFile)
  sampleDF <- read.table(sampleDGEFile,header = T,sep = "\t",stringsAsFactors = F)
  rownames(sampleDF) <- sampleDF[,1]
  sampleDF <- sampleDF[,c(2:ncol(sampleDF))]
  sampleDF$means <- rowMeans(sampleDF)
  sampleDF$genes <- rownames(sampleDF)
  gCol <- ncol(sampleDF)
  mCol <- gCol - 1
  
  newDF <- sampleDF[,c(gCol,mCol)]
  colnames(newDF) <- c("genes",sample)
  
  resultDF <- merge(resultDF,newDF,by=c("genes","genes"),all=TRUE)
}

resultDF <- resultDF[,c(1,3:ncol(resultDF))]
print(head(resultDF))

write.table(x = resultDF,file = outputFile ,sep = "\t")

