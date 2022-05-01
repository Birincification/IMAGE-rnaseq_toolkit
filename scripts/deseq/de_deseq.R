#!/usr/local/bin/Rscript
# call: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <de.method> <out.file>
#
# choose <de.method> out of {'limma', 'edgeR', 'DESeq'}

#print(commandArgs())
#if(length(commandArgs()) != 10) message("usage: Rscript de_rseq.R <exprs.file> <pdat.file> <fdat.file> <de.method> <out.file>")
#stopifnot(length(commandArgs()) == 10)

message("Loading EnrichmentBrowser")
print(.libPaths())
suppressWarnings(suppressPackageStartupMessages(library(EnrichmentBrowser)))

args <- commandArgs(trailingOnly=TRUE)
exprs.file <- args[1]
pdat.file <- args[2]
fdat.file <- args[3]
de.method <- args[4]
out.file <- args[5]
print(out.file)
print(exprs.file)


message("Reading data ...")
eset <- readSE(exprs.file, pdat.file, fdat.file)

message("DE analysis ...")
eset <- deAna(eset, de.method=de.method, padj.method="none")

de.tbl <- rowData(eset)
de.tbl <- de.tbl[, c("ENTREZID", "FC", "PVAL")]
colnames(de.tbl) <- c("GENE.ID", "log2FC", "RAW.PVAL")
de.tbl$ADJ.PVAL <- p.adjust(de.tbl$RAW.PVAL, method="BH")

write.table(de.tbl, file=out.file, row.names=FALSE, quote=FALSE, sep="\t")
message("DE table written to ", out.file)
