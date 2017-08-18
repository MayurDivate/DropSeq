# cell selction plot
# This script will help to plot drop seq data
# in order to determine number of cells 
# As this script can not predict number of cell, 
# I suggest you to plot data using different x axis limits 
# to find correct number of cells
# Author@MayurDivate

# *_cell_readcounts.txt.gz file generated in step 1
readcountFile <- "SRR1873277_cell_readcounts.txt.gz"
a <- read.table(readcountFile,header = F,stringsAsFactors = F) 

imgFile <- gsub("_cell_readcounts.txt.gz","_cell_selction.jpg",readcountFile)
mtitle <- gsub("_cell_readcounts.txt.gz","",readcountFile)

a <- read.table(readcountFile,header = F,stringsAsFactors = F) 
x <- cumsum(a$V1)
x <- x / max(x)

# default x axis limit is 1500
# You can consider changing x axis limit to find correct number of cells ( knee point ) 
x_axis_limit <- 1500

jpeg(imgFile,width = 1000, height = 1000,quality = 100,pointsize = 12,res=200)
par(mar=c(4,4,4,1))
plot(1:length(x),x,type="l",col="blue",xlab = "cell barcodes sorted by number of reads [descending]",
     ylab="cumulative fraction of reads", main = mtitle, xlim=c(0,x_axis_limit))

# change nc to appropiate cell number 
nc <- 570
pointBy <- x[nc]
pointBy
lines(x=c(nc,nc),y=c(min(x),pointBy),col="red")

nctext = paste("cells = ",nc,sep = "")
mtext(nctext,side = 3,adj = 0.9,line = -2)
dev.off()


# done

