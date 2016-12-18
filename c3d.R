#!/usr/bin/Rscript

# Copyright 2016 Mathieu Lupien

# This file is part of C3D.

# C3D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# C3D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with C3D.  If not, see <http://www.gnu.org/licenses/>.

# Cross Cell-type Correlation in DNaseI hypersensitivity (C3D)
# Written by: Tahmid Mehdi 
# Princess Margaret Cancer Centre - University Health Network, December 18, 2016
# Tested on R 3.3.1

# Takes correlations between open regions of chromatin based on DNaseI hypersensitivity signals
# Regions with high correlations are candidates for 3D interactions
# Performs association tests on each candidate & adjusts p-values
# Produces interaction landscapes and tracks in PDF format

args <- commandArgs(trailingOnly = TRUE)
refMapDir <- args[1] # directory of mapped files
outDir <- sub('/$', '', args[2]) # output directory
anchor <- args[3] # anchor file (bed)
bg <- args[4] # file with list of bg files 
window <- args[5] # flanking bps from anchor to search for open regions
rcut <- as.numeric(args[6]) # correlation threshold
pcut <- as.numeric(args[7]) # p-value threshold
qcut <- as.numeric(args[8]) # q-value threshold
corMethod <- tolower(args[9]) # correlation coefficient (pearson, spearman, kendall)
signalMatrixFile <- args[10] # signal data (if available)
figures <- tolower(args[11]) # 'y' if you want interaction landscapes outputted
figureWidth <- args[12] # flanking bps from anchor to display on figure
zoom <- args[13] # how many flanking bps for zoomed landscapes
colour <- args[14] # colours for the interaction plots
tracks <- args[15] # 'y' if you want tracks outputted
sampleName <- args[16]
trackNumber <- as.numeric(args[17])
numSamples <- as.numeric(args[18])
assembly <- args[19] # reference genome
date <- args[20] # timestamp for output pdf
workingDir <- args[21] # working directory

setwd(workingDir)
suppressMessages(library("GenomicRanges"))
suppressMessages(library("Sushi"))
suppressMessages(library("data.table"))
suppressMessages(library("preprocessCore"))
suppressMessages(library("dynamicTreeCut"))

# PRE: str
# POST: GRange
# converts bed file to a granges object
bed_to_granges = function(file) {
  df <- read.table(file, header=F, stringsAsFactors=F)[,1:3]
  names(df) <- c('chr','start','end')
  gr <- with(df, GRanges(chr, IRanges(start, end)))
  return(gr)
}

# PRE: int[>=1] int[>=1] matrix(num)
# POST: vector(num)
# returns correlation and p-value of xth & yth rows of A
getCor = function(x, y, A) {
  if (sd(A[x,])==0 || sd(A[y,])==0) { # avoid div by 0
    return(c(0,NA))
  }
  corr <- cor.test(A[x,], A[y,], method=corMethod, alternative="two.sided")
  return(c(unname(corr$estimate), corr$p.value))
}

# PRE: df(chr1,start1,end1,chr2,start2,end2,score,q_Value) str int[>=1] int[>=1]
# POST: 0<=num<=1
# returns max score of int on chromosome chr between start & end
get_max_cor = function(int, chr, start, end) {
  filteredInt <- int[with(int, chrom1==chr & end1>=start & start1<=end &
                            chrom2==chr & end2>=start & start2<=end),]
  return(max(filteredInt$score))
}

# PRE: str int[>=1]
# POST: vect
# Breaks a string into a vector of strings separated at commas
commaSepStr_to_vector = function(commaSepStr, desiredLength) {
  vect <- unlist(strsplit(commaSepStr, split=","))
  # if the string didn't have desiredLength items,
  # then just make a vector with the first element repeated desiredLength times
  if (length(vect)!=desiredLength) {
    vect <- rep(vect[1], desiredLength)
  }
  return(vect)
}

# PRE: str
# POST: int or str
# converts strands into numbers for the Sushi package
strand_to_num = function(strand) {
  if (strand=="+") {
    return(1)
  } else if (strand=="-") {
    return(-1)
  } else {
    return(".")
  }
}

# PRE: int[>=0] vect[hex]
# POST: str
# Converts the xth element of pal (a palette) to RGB values
getRGB <- function(x, pal) {
  rgbValues <- col2rgb(pal[x])
  return(paste(rgbValues[1],rgbValues[2],rgbValues[3], sep=","))
}

# PRE: matrix or vector
# POST: vector
# Calculates the mean of each row
rowAvg = function(A) {
  if (is.matrix(A)) {
    return(rowMeans(A))
  } else {
    return(A)
  }
}

# if signalMatrixFile not provided, generate matrix from mapped files
if (signalMatrixFile=="") {
  # list of files in refMapDir ending in .map.bed
  refMapDir <-  sub('/$', '', refMapDir)
  refMapFiles <- list.files(refMapDir, pattern="*.map.bed", full.names=TRUE)
  # granges object with regions from reference file
  ref.bed <- bed_to_granges(refMapFiles[1])
  # format names of regions from reference file
  regionNames <- paste(as.character(seqnames(ref.bed)),":",start(ranges(ref.bed)),
                       "-",end(ranges(ref.bed)), sep="")
  
  cat("Merging Mapped Files...\n")
  # merge scores from bg files
  signals <- do.call(cbind, lapply(refMapFiles, 
                                   function(f) read.table(f,header=FALSE, sep="\t")[,4]))
  rownames(signals) <- regionNames
  write.table(signals, paste(outDir,"signalMatrix.txt", sep="/"), 
              row.names=T, col.names=F, quote=F) # write matrix to file
} else { # if signalMatrixFile provided, use that as the matrix
  # load matrix
  signals <- as.matrix(read.table(signalMatrixFile, header=FALSE, row.names=1))
  regionNames <- rownames(signals)
  # create granges object from matrix file
  partsOfBed <- strsplit(regionNames, ":|-")
  chr <- sapply(partsOfBed, function(x) x[1])
  start <- sapply(partsOfBed, function(x) as.numeric(x[2]))
  end <- sapply(partsOfBed, function(x) as.numeric(x[3]))
  df <- data.frame(chr, start, end)
  ref.bed <- with(df, GRanges(chr, IRanges(start, end)))
}

if (! is.numeric(signals)) { # if matrix has NAs
  stop("There are null values in the signal matrix. Check mapped files.\n")
}

signals.norm <- normalize.quantiles(signals) # quantile normalize the signals
rownames(signals.norm) <- regionNames
# calculate distances between samples (1-correlation)
crossGenomeCors <- as.dist(1-cor(signals.norm, method="pearson"))
dendro <- hclust(crossGenomeCors) # hierarchically cluster the samples
# run the dynamic tree cut algorithm to detect clusters
clusters <- as.character(cutreeDynamic(dendro, minClusterSize=1, verbose=0, 
                                       distM=as.matrix(crossGenomeCors), deepSplit=4))
# average normalized signal at each open region for each cluster
avg.norm.signals <- do.call(cbind, lapply(unique(clusters), 
                                          function(p) rowAvg(signals.norm[,which(clusters==p)])))
cat("Clustered samples into",length(unique(clusters)),"clusters\n" )
# vector of IDs in same order as anchor file
ids <- read.table(anchor, header=F, stringsAsFactors=F)[,5]
anchorStrands <- read.table(anchor, header=F, stringsAsFactors=F)[,4]
# convert anchor to GRanges object
anchor.bed <- bed_to_granges(anchor)

# find promoters which overlap reference regions
promoterOverlaps <- findOverlaps(anchor.bed, ref.bed)
if (window != "genome") {
  window <- as.numeric(window)
  # find reference regions within window of each promoter
  leftOverlaps <- findOverlaps(flank(anchor.bed, window), ref.bed)
  rightOverlaps <- findOverlaps(flank(anchor.bed, window, start=FALSE), ref.bed)
}

# initialize list of interaction candidate
anchorStats <- vector(mode="list", length=length(anchor.bed))
file <- paste(outDir,"/results_",date,".txt", sep="")
# change the first letter in corMethod to uppercase for the arc diagrams
corMethodUpper <- paste(toupper(substr(corMethod, 1, 1)), substr(corMethod, 2, nchar(corMethod)), sep="")
resultsHeader <- paste("COORD_1\tCOORD_2\tR_",corMethodUpper,"\tID\tp_Value\tq_Value" , sep="")
junk <- cat(resultsHeader, file=file, append=FALSE, sep="\n")

for (r in 1:length(anchor.bed)) { # iterate through promoters
  COORD_1 <- c() # regions in promoter
  COORD_2 <- c() # distal regions from promoter
  Correlation <- c() # correlation
  p_Value <- c() # p-value
  q_Value <- c() # q-value
  ID <- c() # ID of the anchors
  # regions in rth promoter
  regionsInPromoters <- subjectHits(promoterOverlaps[which(queryHits(promoterOverlaps)==r)])
  if (window=="genome") { # search the entire genome
    candidateRegions <- setdiff(1:length(ref.bed), regionsInPromoters)
  } else {
    regionsInLeft <- subjectHits(leftOverlaps[which(queryHits(leftOverlaps)==r)])
    regionsInRight <- subjectHits(rightOverlaps[which(queryHits(rightOverlaps)==r)])
    # regions in window
    candidateRegions <- unique(c(regionsInLeft, regionsInRight))
    candidateRegions <- candidateRegions[! candidateRegions %in% regionsInPromoters]
  }
  if (length(regionsInPromoters)==0 || length(candidateRegions)==0) {
    cat('\rCalculating Correlations: Processed anchor', r, 'of', length(anchor.bed))
    next
  }
  # calculate correlations between open regions in promoter & distal regions
  regionIndices <- cbind(rep(regionsInPromoters,each=length(candidateRegions)), 
                         rep(candidateRegions, length(regionsInPromoters)))
  corStats <- apply(regionIndices, 1, function(x) getCor(x[1], x[2], signals))
  # record regionNames of anchor DHSs 
  COORD_1 <- regionNames[regionIndices[,1]]
  # record regionNames of distal DHSs 
  COORD_2 <- regionNames[regionIndices[,2]]
  # record correlations, p-values & q-values
  Correlation <- corStats[1,]
  p_Value <- corStats[2,]
  ID <- rep(ids[r], nrow(regionIndices)) # populate ID
  q_Value <- p.adjust(p_Value, method="BH")
  interactionCandidates <- data.frame(COORD_1, COORD_2, Correlation, ID, p_Value, q_Value)
  if (nrow(interactionCandidates)>0) { # stop if genomeStats is empty
    junk <- write.table(interactionCandidates, file=file, sep="\t", col.names=F,row.names=F,quote=F, append=TRUE)
  }
  anchorStats[[r]] <- interactionCandidates
  cat('\rCalculating Correlations: Processed anchor', r, 'of', length(anchor.bed))
}
# concatenate anchorStats data.frames
genomeStats <- rbindlist(anchorStats)
# filter genomeStats
genomeStats <- genomeStats[with(genomeStats, Correlation>=rcut & p_Value<=pcut & q_Value<=qcut),]
cat("\n")
# get ref.bed indices corresponding to COORD_1 & COORD_2
regionIndices.Coord1 <- sapply(genomeStats$COORD_1,
                               function(x) which(regionNames==x))
regionIndices.Coord2 <- sapply(genomeStats$COORD_2,
                               function(x) which(regionNames==x))
if (length(regionIndices.Coord1)<=0 && figures=="y") {
  cat("No interaction candidates passed the filters; No figures generated\n")
}
# get widths for figures & tracks
figureWidth <- as.numeric(commaSepStr_to_vector(figureWidth, length(anchor.bed)))
# extract zoom lengths
zoom <- as.numeric(commaSepStr_to_vector(zoom, length(anchor.bed)))
# graphics --------------------------------------------------------------------
if (figures=="y" && length(regionIndices.Coord1)>0) {
  # create data frame for anchors
  chr <- as.character(seqnames(anchor.bed))
  start <- start(ranges(anchor.bed))
  end <- end(ranges(anchor.bed))
  name <- sapply(strsplit(ids, "_"), function(x) head(x, n=1))
  strand <- sapply(anchorStrands, function(x) strand_to_num(x)) 
  score <- rep(".", length(name))
  p_Value <- q_Value <- rep(0, length(name))
  
  anchor.df <- data.frame(chr,start,end, strand, score, p_Value, q_Value, name, row=1, color="purple")
  
  # extract colours for interactions
  colour <- commaSepStr_to_vector(colour, 4)
  # make interaction data frame 
  chrom1 <- as.character(seqnames(ref.bed[regionIndices.Coord1]))
  start1 <- start(ranges(ref.bed[regionIndices.Coord1]))
  end1 <- end(ranges(ref.bed[regionIndices.Coord1]))
  chrom2 <- as.character(seqnames(ref.bed[regionIndices.Coord2]))
  start2 <- start(ranges(ref.bed[regionIndices.Coord2]))
  end2 <- end(ranges(ref.bed[regionIndices.Coord2]))
  color <- rep("black", length(genomeStats$q_Value))
  # make colours based on q-values
  color[which(genomeStats$q_Value>0.05)] <- colour[1]
  color[which(0.01<genomeStats$q_Value & genomeStats$q_Value<=0.05)] <- colour[2]
  color[which(0.001<genomeStats$q_Value & genomeStats$q_Value<=0.01)] <- colour[3]
  color[which(genomeStats$q_Value<=0.001)] <- colour[4]
  interactions.bedpe <- data.frame(chrom1, start1, end1, chrom2, start2, end2, 
                                   score=genomeStats$Correlation, q_Value=genomeStats$q_Value, color)
  
  zoomOverWidth <- which(zoom>figureWidth)
  zoom[zoomOverWidth] <- figureWidth[zoomOverWidth] # truncate values over window to window
  noZoom <- which(zoom<=0) # figures with no zoom
  
  corName <- paste(toupper(substr(corMethod, 1, 1)), 
                   substr(corMethod, 2, nchar(corMethod)), sep="")
  # make pdf file with figures
  pdf(file=paste(outDir,"/figures_",date,".pdf", sep=""), height=11, width=8.5)
  for (r in 1:length(anchor.bed)) { # make figure for each anchor
    # set up the region
    chrom = as.character(seqnames(anchor.bed[r]))
    chromstart <- start(ranges(anchor.bed[r])) - figureWidth[r]
    chromend <- end(ranges(anchor.bed[r])) + figureWidth[r]
    # turn off warnings temporarily
    defaultWarn <- getOption("warn")
    options(warn = -1)
    maxCor <- get_max_cor(interactions.bedpe, chrom,chromstart,chromend)
    options(warn = defaultWarn)
    if (r %in% noZoom) { # figure with no zoom
      layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2), 7,2, byrow=TRUE))
      par(mar=c(4,6,1,8), oma=rep(0.75,4), xpd=TRUE) 
      # landscape plot
      intLandscape = plotBedpe(interactions.bedpe, chrom,chromstart,chromend, 
                               heights=interactions.bedpe$score, 
                               lwdby=0.5*dnorm(interactions.bedpe$q_Value, mean=0, sd=0.1), lwdrange=c(0,2),
                               plottype="loops", color=interactions.bedpe$color, ymax=1/maxCor)
      labelgenome(chrom,chromstart,chromend, n=5,scale="bp", 
                  chromline=0.25,scaleline=0.25)
      legend("right",inset=-0.16,title=expression(bold("q-value")),
             legend=c(expression(""<="1"),expression(""<="0.05"),
                      expression(""<="0.01"),expression(""<="0.001")), 
             col=colour,lty=1,lwd=2.5,text.font=2,cex=1.25)
      axis(side=2,las=2,at=seq(0,1,0.1), xpd=TRUE)
      mtext(paste(corName,"Correlation",sep=" "),side=2,line=2.5,font=2)
      # tracks for anchors
      par(mar=c(1,6,1,8)) 
      anchorTrack = plotBed(beddata=anchor.df,chrom=chrom,chromstart=chromstart,chromend=chromend,
                            rownumber=anchor.df$row, type="region", color=anchor.df$color,row="given",
                            rowlabels=c("Anchors"), rowlabelcol="black", rowlabelcex=1.25)
    } else { # figure with zoom
      layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3), 7,2, byrow=TRUE))
      par(mar=c(4,6,1,8), oma=rep(0.75,4), xpd=TRUE) 
      # landscape plot
      intLandscape = plotBedpe(interactions.bedpe, chrom,chromstart,chromend, 
                               heights=interactions.bedpe$score, 
                               lwdby=0.5*dnorm(interactions.bedpe$q_Value, mean=0, sd=0.1), lwdrange=c(0,2),
                               plottype="loops", color=interactions.bedpe$color, ymax=1/maxCor)
      labelgenome(chrom,chromstart,chromend, n=5,scale="bp", 
                  chromline=0.25,scaleline=0.25)
      legend("right",inset=-0.16,title=expression(bold("q-value")),
             legend=c(expression(""<="1"),expression(""<="0.05"),
                      expression(""<="0.01"),expression(""<="0.001")), 
             col=colour,lty=1,lwd=2.5,text.font=2,cex=1.25)
      axis(side=2,las=2,at=seq(0,1,0.1), xpd=TRUE)
      mtext(paste(corName,"Correlation",sep=" "),side=2,line=2.5,font=2)
      # zoomed region
      par(mar=c(1,6,1,8)) 
      regionZoom <- c(start(ranges(anchor.bed[r]))-zoom[r], end(ranges(anchor.bed[r]))+zoom[r])
      zoomsregion(regionZoom,extend=c(-1,0.13),wideextend=0,offsets=c(0,0))
      options(warn = -1)
      maxCor <- get_max_cor(interactions.bedpe, chrom,regionZoom[1],regionZoom[2])
      options(warn = defaultWarn)
      # zoomed landscape plot
      intLandscapeZoom = plotBedpe(interactions.bedpe, chrom,chromstart=regionZoom[1],chromend=regionZoom[2], 
                                   heights=interactions.bedpe$score, 
                                   lwdby=0.5*dnorm(interactions.bedpe$q_Value, mean=0, sd=0.1), lwdrange=c(0,2),
                                   plottype="loops", color=interactions.bedpe$color, ymax=1/maxCor)
      zoombox()
      labelgenome(chrom,chromstart=regionZoom[1],chromend=regionZoom[2], n=5,scale="bp", 
                  chromline=0.25,scaleline=0.25)
      axis(side=2,las=2,at=seq(0,1,0.1), xpd=TRUE)
      mtext(paste(corName,"Correlation",sep=" "),side=2,line=2.5,font=2)
      # tracks for anchors
      par(mar=c(1,6,1,8)) 
      anchorTrack = plotBed(beddata=anchor.df,chrom=chrom,chromstart=regionZoom[1],chromend=regionZoom[2],
                            rownumber=anchor.df$row, type="region", color=anchor.df$color,row="given",
                            rowlabels=c("Anchors"), rowlabelcol="black", rowlabelcex=1.25)
    }
    mtext(paste("Figure ",r,": Interaction Landscape of", ids[r], sep=""), side=1)
    cat('\rGenerating Figures: Created figure', r, 'of', length(anchor.bed))
  }
  cat("\n")
  turnDevOff <- dev.off()
}

# tracks ----------------------------------------------------------------------
if (tracks=="y") {
  trackPal <- colorRampPalette(c("blue", "purple", "red", "orange"))(numSamples)
  for (r in 1:length(anchor.bed)) {
    if (is.null(anchorStats[[r]])) {
      cat('\rGenerating Tracks: Created track', r, 'of', length(anchor.bed))
      next
    }
    file <- paste(outDir,"/",ids[r],".anchor", sep="")
    # create & print the browser header to a file
    browserHeader <- paste("browser position ",as.character(seqnames(anchor.bed[r])),":",
                           start(anchor.bed[r])-figureWidth[r],"-",end(anchor.bed[r])+figureWidth[r],
                           "\nbrowser hide all\nbrowser pack refGene\nbrowser full altGraph", sep="")
    junk <- cat(browserHeader, file=file, append=FALSE, sep="\n")
    # create & print the track header for the anchor
    trackHeader <- paste("track type=bedGraph name=\"Anchor\" description=\"",gsub("_.*", "", ids[r]),
                         "\" visibility=full color=0,0,0 altColor=0,0,0 viewLimits=0:1 autoScale=off gridDefault=on db=", assembly, sep="")
    junk <- cat(trackHeader, file=file, append=TRUE, sep="\n")
    # write the anchor track
    anchorLine <- paste(as.character(seqnames(anchor.bed[r])),start(anchor.bed[r]),end(anchor.bed[r]),1, sep="\t")
    junk <- cat(anchorLine, file=file, append=TRUE, sep="\n")
    
    file <- paste(outDir,"/",ids[r],".bedGraph", sep="")
    # create & print the track header for the distal DHSs
    trackHeader <- paste("track type=bedGraph name=\"",sampleName,"\" description=\" \" visibility=full color=",getRGB(trackNumber, trackPal), 
                         " altColor=0,0,0 viewLimits=0:1 autoScale=off gridDefault=on yLineMark=",rcut," yLineOnOff=on db=", assembly, sep="")
    junk <- cat(trackHeader, file=file, append=FALSE, sep="\n")
    # create a temporary data frame for storing distal DHSs & their correlations with the anchor
    tmp.df <- data.frame(distal=anchorStats[[r]]$COORD_2, corr=anchorStats[[r]]$Correlation)
    # if a distal DHS appears more than once, pick the one with the highest correlation
    # this will only happen if the promoter has multiple open regions
    tmp.df <- tmp.df[with(tmp.df, order(distal, -corr)), ]
    tmp.df <- tmp.df[ !duplicated(tmp.df$distal), ]
    tmp.df <- tmp.df[tmp.df$corr>=0, ] # filter out negative correlations
    # split the distal DHSs by chr, start & stop
    distalDHS <- strsplit(as.character(tmp.df$distal), ":|-")
    chr <- sapply(distalDHS, function(x) x[1])
    start <- sapply(distalDHS, function(x) x[2])
    end <- sapply(distalDHS, function(x) x[3])
    correlation <- round(tmp.df$corr,2) # round correlations to 2 decimals
    # write track to file
    track <- data.frame(chr, start, end, correlation)
    junk <- write.table(track, file=file, sep="\t", col.names=F,row.names=F,quote=F, append=TRUE)
    cat('\rGenerating Tracks: Created track', r, 'of', length(anchor.bed))
  }
  cat("\n")
}
cat("Done\n")
