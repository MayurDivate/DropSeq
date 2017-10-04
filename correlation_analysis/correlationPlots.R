# this script will give the matrix of Avergae genes per million
# for all the samples in specified folder 
# author@MayurDivate
# 


# specifiy the folder path
setwd("~/Work/0_GitHub/DropSeq/correlation_analysis/sample_data/")

# file containing matrix of avergae transcript per gene
fileAT <- "AverageTrancriptsPerGene.txt"

# file containing total aligned reads per sample 
fileTotalReads <- "TotalReadsPersample.txt"

humanAT <- read.table(fileAT,header = T,sep = "\t",stringsAsFactors = F)
humanAT[is.na(humanAT)] = 0
head(humanAT[is.na(humanAT[1,]),])

row.names(humanAT) <- humanAT[,1]
humanAT <- humanAT[,c(2:ncol(humanAT))]

# plot every sample with last sample
nsamples <- ncol(humanAT) - 1
samples <- colnames(humanAT)

totalReaddf <- read.table(fileTotalReads,header=F,stringsAsFactors = F,sep = "\t")
head(totalReaddf)
head(humanAT)

for(i in 1:nsamples){
  df <- humanAT[,c(i,ncol(humanAT))]
    dfX <- df[,c(1,nsamples+1)]
    head(dfX)
    dfX <- dfX[rowSums(dfX) > 0,]
    options(digits = 12)
    dfX <- round(dfX[,c(1:2)],4)

    scaleA <- totalReaddf[i,2]
    scaleA <-  scaleA / 1000000
    scaleA <- round(scaleA,2) 
    scaleA
    
    scaleB <- totalReaddf[nsamples+(x-1),2]
    scaleB <- scaleB / 1000000
    scaleB <- round(scaleB,2) 
    scaleB
    
    head(dfX)
    dfX[,1] <- dfX[,1]*scaleA
    head(dfX)
    dfX[,2] <- dfX[,2]*scaleB
    head(dfX)
    dfX <- round(dfX[,c(1:2)],4)
    dfX <- dfX+1
    head(dfX)
    
    corText <- cor( dfX[,1] , dfX[,2], method = "pearson" )
    corText <- round(corText,2)
    corText <- paste("pearson R=",corText)
    corText 
    
    imgFile <- paste(samples[i], samples[nsamples+(x-1)],sep = "_vs_");
    imgFile <- paste(imgFile,"_pearson.jpg",sep = "")
    print(imgFile)
    print(i)
    print(x)
    
    
    dfX <- dfX + 1
    jpeg(imgFile,width = 1000, height = 1000,quality = 100,pointsize = 12,res=200)
    par(mar=c(5,5,2,1))
    plot(log10(dfX[,1]),log10(dfX[,2]),pch=16, xlab = samples[i], ylab = samples[nsamples+1])
    mtext(corText)  
    abline(lm(log10(dfX[,2])~log10(dfX[,1])))
    dev.off()
    
    plot(log10(dfX[,1]),log10(dfX[,2]),pch=16, xlab = samples[i], ylab = samples[nsamples+(x-1)])
    plot <- ggplot(dfX,aes(x = log10(dfX[,1]), y= log10(dfX[,2]))) + geom_point()
    plot    
    
}


