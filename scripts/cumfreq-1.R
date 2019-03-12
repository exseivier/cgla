#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)
input <- args[1]
output <- gsub(".cat$", ".pdf", input)
#rdata <- args[2] this is for version 2
#rd_obj <- args[3]

load(".RData")
#table <- rd_obj$table

GSNO_YPD.motif <- GSNO_YPD.table
kmer <- read.delim(input, header=T, row.names=1)
rownames(kmer) <- gsub(">", "", rownames(kmer))
kmer <- kmer[rownames(kmer) %in% rownames(GSNO_YPD.motif),]
GSNO_YPD.motif <- GSNO_YPD.motif[rownames(GSNO_YPD.motif) %in% rownames(kmer),]

if(length(rownames(kmer)) == length(rownames(GSNO_YPD.motif))){
	print("Length kmers and DGE data match..")
}


kmer <- kmer[order(rownames(kmer)),]
GSNO_YPD.motif <- GSNO_YPD.motif[order(rownames(GSNO_YPD.motif)),]
GSNO_YPD.motif <- cbind(GSNO_YPD.motif, kmer)
GSNO_YPD.motif <- GSNO_YPD.motif[order(-GSNO_YPD.motif$logFC),]
GSNO_YPD.motif.vn <- GSNO_YPD.motif[GSNO_YPD.motif$very_near > 0,]
GSNO_YPD.motif.n <- GSNO_YPD.motif[GSNO_YPD.motif$near > 0,]
GSNO_YPD.motif.m <- GSNO_YPD.motif[GSNO_YPD.motif$moderate > 0,]
GSNO_YPD.motif.f <- GSNO_YPD.motif[GSNO_YPD.motif$far > 0,]
GSNO_YPD.motif.vf <- GSNO_YPD.motif[GSNO_YPD.motif$very_far > 0,]


GSNO_YPD.motif.brk <- seq(min(GSNO_YPD.motif$logFC), max(GSNO_YPD.motif$logFC), by = 0.5)
GSNO_YPD.motif.vn.brk <- seq(min(GSNO_YPD.motif.vn$logFC), max(GSNO_YPD.motif.vn$logFC), by = 0.5)
GSNO_YPD.motif.n.brk <- seq(min(GSNO_YPD.motif.n$logFC), max(GSNO_YPD.motif.n$logFC), by = 0.5)
GSNO_YPD.motif.m.brk <- seq(min(GSNO_YPD.motif.m$logFC), max(GSNO_YPD.motif.m$logFC), by = 0.5)
GSNO_YPD.motif.f.brk <- seq(min(GSNO_YPD.motif.f$logFC), max(GSNO_YPD.motif.f$logFC), by = 0.5)
GSNO_YPD.motif.vf.brk <- seq(min(GSNO_YPD.motif.vf$logFC), max(GSNO_YPD.motif.vf$logFC), by = 0.5)

GSNO_YPD.motif.cut <- cut(GSNO_YPD.motif$logFC, GSNO_YPD.motif.brk)
GSNO_YPD.motif.vn.cut <- cut(GSNO_YPD.motif.vn$logFC, GSNO_YPD.motif.vn.brk)
GSNO_YPD.motif.n.cut <- cut(GSNO_YPD.motif.n$logFC, GSNO_YPD.motif.n.brk)
GSNO_YPD.motif.m.cut <- cut(GSNO_YPD.motif.m$logFC, GSNO_YPD.motif.m.brk)
GSNO_YPD.motif.f.cut <- cut(GSNO_YPD.motif.f$logFC, GSNO_YPD.motif.f.brk)
GSNO_YPD.motif.vf.cut <- cut(GSNO_YPD.motif.vf$logFC, GSNO_YPD.motif.vf.brk)

GSNO_YPD.motif.freq <- table(GSNO_YPD.motif.cut)
GSNO_YPD.motif.vn.freq <- table(GSNO_YPD.motif.vn.cut)
GSNO_YPD.motif.n.freq <- table(GSNO_YPD.motif.n.cut)
GSNO_YPD.motif.m.freq <- table(GSNO_YPD.motif.m.cut)
GSNO_YPD.motif.f.freq <- table(GSNO_YPD.motif.f.cut)
GSNO_YPD.motif.vf.freq <- table(GSNO_YPD.motif.vf.cut)

GSNO_YPD.motif.cumfreq <- c(0, cumsum(GSNO_YPD.motif.freq))
GSNO_YPD.motif.vn.cumfreq <- c(0, cumsum(GSNO_YPD.motif.vn.freq))
GSNO_YPD.motif.n.cumfreq <- c(0, cumsum(GSNO_YPD.motif.n.freq))
GSNO_YPD.motif.m.cumfreq <- c(0, cumsum(GSNO_YPD.motif.m.freq))
GSNO_YPD.motif.f.cumfreq <- c(0, cumsum(GSNO_YPD.motif.f.freq))
GSNO_YPD.motif.vf.cumfreq <- c(0, cumsum(GSNO_YPD.motif.vf.freq))

GSNO_YPD.motif.vn.cumfreq <- GSNO_YPD.motif.vn.cumfreq/max(GSNO_YPD.motif.vn.cumfreq)
GSNO_YPD.motif.n.cumfreq <- GSNO_YPD.motif.n.cumfreq/max(GSNO_YPD.motif.n.cumfreq)
GSNO_YPD.motif.m.cumfreq <- GSNO_YPD.motif.m.cumfreq/max(GSNO_YPD.motif.m.cumfreq)
GSNO_YPD.motif.f.cumfreq <- GSNO_YPD.motif.f.cumfreq/max(GSNO_YPD.motif.f.cumfreq)
GSNO_YPD.motif.vf.cumfreq <- GSNO_YPD.motif.vf.cumfreq/max(GSNO_YPD.motif.vf.cumfreq)
GSNO_YPD.motif.cumfreq <- GSNO_YPD.motif.cumfreq/max(GSNO_YPD.motif.cumfreq)

lo <- loess(GSNO_YPD.motif.cumfreq ~ GSNO_YPD.motif.brk)
lo.vn <- loess(GSNO_YPD.motif.vn.cumfreq ~ GSNO_YPD.motif.vn.brk)
lo.n <- loess(GSNO_YPD.motif.n.cumfreq ~ GSNO_YPD.motif.n.brk)
lo.m <- loess(GSNO_YPD.motif.m.cumfreq ~ GSNO_YPD.motif.m.brk)
lo.f <- loess(GSNO_YPD.motif.f.cumfreq ~ GSNO_YPD.motif.f.brk)
lo.vf <- loess(GSNO_YPD.motif.vf.cumfreq ~ GSNO_YPD.motif.vf.brk)

lfc <- sample(GSNO_YPD.motif$logFC, size=length(GSNO_YPD.motif.vn$logFC))
w.vn <- round(wilcox.test(lfc, GSNO_YPD.motif.vn$logFC)$p.value, 4)
lfc <- sample(GSNO_YPD.motif$logFC, size=length(GSNO_YPD.motif.n$logFC))
w.n <- round(wilcox.test(GSNO_YPD.motif$logFC, GSNO_YPD.motif.n$logFC)$p.value, 4)
lfc <- sample(GSNO_YPD.motif$logFC, size=length(GSNO_YPD.motif.m$logFC))
w.m <- round(wilcox.test(GSNO_YPD.motif$logFC, GSNO_YPD.motif.m$logFC)$p.value, 4)
lfc <- sample(GSNO_YPD.motif$logFC, size=length(GSNO_YPD.motif.f$logFC))
w.f <- round(wilcox.test(GSNO_YPD.motif$logFC, GSNO_YPD.motif.f$logFC)$p.value, 4)
lfc <- sample(GSNO_YPD.motif$logFC, size=length(GSNO_YPD.motif.vf$logFC))
w.vf <- round(wilcox.test(GSNO_YPD.motif$logFC, GSNO_YPD.motif.vf$logFC)$p.value, 4)



pdf(output, width=7, height=7)
plot(lo$x, lo$y, type="l", col="black", ylab="Frequency", xlab="logFC", main=gsub(".pdf$", "", output))
lines(lo.vn$x, lo.vn$y, col="red")
lines(lo.n$x, lo.n$y, col="darkred")
lines(lo.m$x, lo.m$y, col="blue")
lines(lo.f$x, lo.f$y, col="darkblue")
lines(lo.vf$x, lo.vf$y, col="darkgreen")
legend("topleft", legend=c("total", paste("very_near ", w.vn, sep=""), paste("near ", w.n, sep=""), paste("moderate ", w.m, sep=""), paste("far ", w.f, sep=""), paste("very_far ", w.vf, sep="")), col=c("black", "red", "darkred", "blue", "darkblue", "darkgreen"), lwd=2, cex=0.7)
dev.off()
