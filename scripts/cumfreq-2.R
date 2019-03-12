#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)
input <- args[1] # *.cat file: kmer matching bins along the sequence. This file is the output of mathcingPattern-2.py
print(input)
output <- gsub(".cat$", ".pdf", input)
print(output)
rdata <- args[2] # RData file produced by edgeR contrast analysis
print(rdata)
rd_obj <- args[3] # The name of the contrast analysis DGE_fit edgeR object
print(rd_obj)

load(rdata)
rd_obj <- eval(parse(text=rd_obj))
table <- rd_obj$table # Extract the table of the contrast analysis

motif <- table
kmer <- read.delim(input, header=T, row.names=1) # Reading the *.cat file and allocating info to kmer data frame
rownames(kmer) <- gsub(">", "", rownames(kmer))
kmer <- kmer[rownames(kmer) %in% rownames(motif),] # Ensure that kmer rownames are in motif rownames
motif <- motif[rownames(motif) %in% rownames(kmer),] # Esnure that motif rownames are in the filtered kmer rownames

# rownames length for kmer and motif are double checked
if(length(rownames(kmer)) == length(rownames(motif))){
	print("Length kmers and DGE data match..")
}

# ordering kmer by gene name
kmer <- kmer[order(rownames(kmer)),]
# ordering motif by gene name
motif <- motif[order(rownames(motif)),]
# binding kmer and motif by columns
motif <- cbind(motif, kmer)
# ordering motif by logFC from higher to lower values
motif <- motif[order(-motif$logFC),]
# Selecting sub-dataset which members have at least one match in a very near location in promoters
motif.vn <- motif[motif$very_near > 0,]
# Selecting sub-dataset which members have at least one match in a near location in promoters
motif.n <- motif[motif$near > 0,]
# Selecting sub-dataset which members have at least one match in the middle location in promoters
motif.m <- motif[motif$moderate > 0,]
# Selecting sub-dataset which members have at least one match in a far location in promoters
motif.f <- motif[motif$far > 0,]
# Selecting sub-dataset which members have at least one match in a very far location in promoters
motif.vf <- motif[motif$very_far > 0,]

# Preparing the breaks for every dataset
motif.brk <- seq(min(motif$logFC), max(motif$logFC), by = 0.01)
motif.vn.brk <- seq(min(motif.vn$logFC), max(motif.vn$logFC), by = 0.01)
motif.n.brk <- seq(min(motif.n$logFC), max(motif.n$logFC), by = 0.01)
motif.m.brk <- seq(min(motif.m$logFC), max(motif.m$logFC), by = 0.01)
motif.f.brk <- seq(min(motif.f$logFC), max(motif.f$logFC), by = 0.01)
motif.vf.brk <- seq(min(motif.vf$logFC), max(motif.vf$logFC), by = 0.01)

# Allocating the logFC values to the corresponding break
motif.cut <- cut(motif$logFC, motif.brk)
motif.vn.cut <- cut(motif.vn$logFC, motif.vn.brk)
motif.n.cut <- cut(motif.n$logFC, motif.n.brk)
motif.m.cut <- cut(motif.m$logFC, motif.m.brk)
motif.f.cut <- cut(motif.f$logFC, motif.f.brk)
motif.vf.cut <- cut(motif.vf$logFC, motif.vf.brk)

# Counting the frequency of logFC vales for every break (bins)
motif.freq <- table(motif.cut)
motif.vn.freq <- table(motif.vn.cut)
motif.n.freq <- table(motif.n.cut)
motif.m.freq <- table(motif.m.cut)
motif.f.freq <- table(motif.f.cut)
motif.vf.freq <- table(motif.vf.cut)

# Performing the cummulative sum for every break (bins)
motif.cumfreq <- c(0, cumsum(motif.freq))
motif.vn.cumfreq <- c(0, cumsum(motif.vn.freq))
motif.n.cumfreq <- c(0, cumsum(motif.n.freq))
motif.m.cumfreq <- c(0, cumsum(motif.m.freq))
motif.f.cumfreq <- c(0, cumsum(motif.f.freq))
motif.vf.cumfreq <- c(0, cumsum(motif.vf.freq))

# Calculating the coefficient for the cummulative frequency for every bins
motif.vn.cumfreq <- motif.vn.cumfreq/max(motif.vn.cumfreq)
motif.n.cumfreq <- motif.n.cumfreq/max(motif.n.cumfreq)
motif.m.cumfreq <- motif.m.cumfreq/max(motif.m.cumfreq)
motif.f.cumfreq <- motif.f.cumfreq/max(motif.f.cumfreq)
motif.vf.cumfreq <- motif.vf.cumfreq/max(motif.vf.cumfreq)
motif.cumfreq <- motif.cumfreq/max(motif.cumfreq)

# Fitting the data to a polynomial model
lo <- loess(motif.cumfreq ~ motif.brk)
lo.vn <- loess(motif.vn.cumfreq ~ motif.vn.brk)
lo.n <- loess(motif.n.cumfreq ~ motif.n.brk)
lo.m <- loess(motif.m.cumfreq ~ motif.m.brk)
lo.f <- loess(motif.f.cumfreq ~ motif.f.brk)
lo.vf <- loess(motif.vf.cumfreq ~ motif.vf.brk)

# Performing wilcoxon rank-sum test for every dataset
lfc <- sample(motif$logFC, size=length(motif.vn$logFC))
w.vn <- round(wilcox.test(lfc, motif.vn$logFC)$p.value, 4)
lfc <- sample(motif$logFC, size=length(motif.n$logFC))
w.n <- round(wilcox.test(motif$logFC, motif.n$logFC)$p.value, 4)
lfc <- sample(motif$logFC, size=length(motif.m$logFC))
w.m <- round(wilcox.test(motif$logFC, motif.m$logFC)$p.value, 4)
lfc <- sample(motif$logFC, size=length(motif.f$logFC))
w.f <- round(wilcox.test(motif$logFC, motif.f$logFC)$p.value, 4)
lfc <- sample(motif$logFC, size=length(motif.vf$logFC))
w.vf <- round(wilcox.test(motif$logFC, motif.vf$logFC)$p.value, 4)

# Plotting the breaks vs. cummulative frequencies
pdf(output, width=7, height=7)
#plot(lo$x, lo$fitted, type="l", col="black", ylab="Frequency", xlab="logFC", main=gsub(".pdf$", "", output))
#lines(lo.vn$x, lo.vn$fitted, col="red")
#lines(lo.n$x, lo.n$fitted, col="blue")
#lines(lo.m$x, lo.m$fitted, col="green")
#lines(lo.f$x, lo.f$fitted, col="yellow")
#lines(lo.vf$x, lo.vf$fitted, col="orange")

plot(lo$x, lo$y, type="l", lwd=3, col="black", ylab="Frequency", xlab="logFC", main=gsub(".pdf$", "", output))
lines(lo.vn$x, lo.vn$y, col="gray60", lty=2, lwd=3)
lines(lo.n$x, lo.n$y, col="gray30", lty=2, lwd=3)
lines(lo.m$x, lo.m$y, col="black", lty=3, lwd=3)
lines(lo.f$x, lo.f$y, col="gray60", lty=4, lwd=3)
lines(lo.vf$x, lo.vf$y, col="gray30", lty=4, lwd=3)
legend("topleft", legend=c("total", paste("very_near ", w.vn, sep=""), paste("near ", w.n, sep=""), paste("moderate ", w.m, sep=""), paste("far ", w.f, sep=""), paste("very_far ", w.vf, sep="")), col=c("black", "gray60", "gray30", "black", "gray60", "gray30"), lwd=2, lty=c(1, 2, 2, 3, 4, 4), cex=0.7)
dev.off()
