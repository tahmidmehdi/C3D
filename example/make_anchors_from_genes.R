# Written by: Paul Guilhamon, Tahmid Mehdi 
# Princess Margaret Cancer Centre - University Health Network, November 29, 2016
# Tested on R 3.3.1

# Creates a BED file of anchors around transcription start sites (TSSs) of specified genes

args <- commandArgs(trailingOnly = TRUE)

genes<-read.delim(as.character(args[1]),header=F,stringsAsFactors=F) # list of genes from file, 1 gene per line
up<-as.numeric(args[2]) # number of bps upstream from TSS
down<-as.numeric(args[3]) # number of bps downstream from TSS
outFile <- as.character(args[4]) # name of output BED file
# read the gencode file which contains gene annotations
gc<-read.delim("gencodev19.txt",stringsAsFactors=F)
# filter gc for desired genes
geneAnnotations <- gc[which(gc$name2 %in% genes[,1]),]
# create data frame for gene annotations for output
geneAnnotations.bed <-data.frame(Loc=character(nrow(geneAnnotations)),Strand=geneAnnotations$strand,Name=character(nrow(geneAnnotations)), stringsAsFactors=F)
# create bed file
for (i in 1:nrow(geneAnnotations)){
  if (geneAnnotations$strand[i]=="+"){
    geneAnnotations.bed$Loc[i]<-paste(geneAnnotations$chrom[i],geneAnnotations$txStart[i]-up,geneAnnotations$txStart[i]+down,sep="\t")
    geneAnnotations.bed$Name[i]<-paste(geneAnnotations$name2[i],"_",geneAnnotations$X.name[i],sep="")
  }
  else if (geneAnnotations$strand[i]=="-"){
    geneAnnotations.bed$Loc[i]<-paste(geneAnnotations$chrom[i],geneAnnotations$txEnd[i]-down,geneAnnotations$txEnd[i]+up,sep="\t")
    geneAnnotations.bed$Name[i]<-paste(geneAnnotations$name2[i],"_",geneAnnotations$X.name[i],sep="")
  }
}
# write bed file
write.table(geneAnnotations.bed, outFile, sep="\t",col.names=F,row.names=F,quote=F)
