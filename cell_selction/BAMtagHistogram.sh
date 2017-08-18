# This script will help to generate input for 
# cell selection plot
# You should have BAMTagHistogram and samtools in system path
# To install BAMTagHistogram follow the link below 
# http://mccarrolllab.com/dropseq/
# or 
# http://mccarrolllab.com/download/922/
# Download dropseq tools and install
# Author@MayurDivate 

in=$1
out=`basename $1 _clean.bam`
cbam=$out"_Dropseq_csrt.bam"
output=$out"_cell_readcounts.txt.gz"

echo "--------------------------------------------------------"

echo "sort bam file to create coordinate sorted bam file "

# add more CPUs for samtools if required using @ parameter  

samtools sort $1 > $cbam 

echo "bam to histogram program"

BAMTagHistogram \
INPUT=$in \
OUTPUT=$output \
TAG=XC \
READ_QUALITY=10

echo "-------------------DONE-------------------------------"

