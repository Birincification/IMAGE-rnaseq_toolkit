#!/usr/local/bin/Rscript
library(argparse)
library(GenomicFeatures)

parser <- ArgumentParser()

parser$add_argument('--gtf')
parser$add_argument('--outdir')
parser$add_argument('--organism')
parser$add_argument('--taxonomyId')

args <- parser$parse_args(commandArgs(trailingOnly=T))


annotation <- makeTxDbFromGFF(args$gtf, format="gtf", organism=args$organism, taxonomyId=as.numeric(args$taxonomyId))
saveDb(annotation, file.path(args$outdir, "txdb.sqlite"))

tx2gene <- transcriptLengths(annotation, T, T, T)
rownames(tx2gene) <- tx2gene$tx_name

tx2gene$tx_name <- NULL
save(tx2gene, file=file.path(args$outdir, "tx2gene.RData"))
