#!/usr/local/bin/Rscript
library(argparse)
library(tximport)
library(DRIMSeq)
library(BiocParallel)
#library(biomaRt)

parser <- ArgumentParser()

parser$add_argument('--counts')
parser$add_argument('--pdata')
parser$add_argument('--tool', default="kallisto", choices=c("kallisto", "salmon"))
parser$add_argument('--seed', default=1984)
parser$add_argument('--ncores', default=4)
parser$add_argument('--outfile')
parser$add_argument('--tx2gene', required=TRUE)
#parser$add_argument('--mart_dataset', default="hsapiens_gene_ensembl")
#parser$add_argument('--mart_version', default=75)
#parser$add_argument('--mart_grch', required=FALSE)


args <- parser$parse_args(commandArgs(trailingOnly=T))

SEED <- as.numeric(args$seed)
set.seed(SEED)
ncores <- as.numeric(args$ncores)

p.data <- read.csv(args$pdata, sep="\t", stringsAsFactors=F)
p.data$sample_id <- p.data$sample
p.data$group <- factor(p.data$condition)
p.data$sample <- NULL
p.data$condition <- NULL
#print(p.data)

if(args$tool == "kallisto")
{
    files <- file.path(args$counts, p.data$sample_id, "abundance.h5")
    print(files)
    txi <- tximport(files, type="kallisto", txOut=TRUE)

} else
{
    message("only kallisto as input is implemented")
    q()
}

counts <- as.data.frame(txi$counts)

colnames(counts) <- p.data$sample_id
counts$feature_id <- rownames(counts)

print(length(counts$feature_id))
print(length(unique(counts$feature_id)))
print(head(counts))

print("imported files...")

#mart <- biomaRt::useEnsembl(biomart="ensembl", GRCh=37, version=75, dataset="hsapiens_gene_ensembl")
#print("loading mart...")
#print(sprintf("datatset: %s", args$mart_dataset))
#print(sprintf("version: %s", args$mart_version))
#if(!is.null(args$mart_grch))
#{
#    print(sprintf("GRCh: %s", args$mart_grch))
#    mart <- biomaRt::useEnsembl(biomart="ensembl", GRCh=as.numeric(args$mart_grch), version=as.numeric(args$mart_version), dataset=args$mart_dataset)
#} else
#{
#    mart <- biomaRt::useEnsembl(biomart="ensembl", version=as.numeric(args$mart_version), dataset=args$mart_dataset)
#}

load(args$tx2gene)


#print("loaded marts...")
#t2g <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"), mart=mart)
#print("loaded t2g")
#rownames(t2g) <- t2g$ensembl_transcript_id

#counts$gene_id <- t2g[counts$feature_id, "ensembl_gene_id"]
counts$gene_id <- tx2gene[counts$feature_id, "gene_id"]

print("transcripts that were not in annotation:")
print(sum(!(counts$feature_id %in% rownames(tx2gene))))


counts <-  counts[rowSums(counts[, 1:nrow(p.data)]) > 0 , ]

print(head(counts))
print(p.data)

d <- dmDSdata(counts=counts, samples=p.data)
print("loaded data")

d <- dmFilter(d, min_samps_gene_expr=nrow(p.data), min_samps_feature_expr=3, min_gene_expr=10, min_feature_expr=10)
print("filtered data")

design_full <- model.matrix(~ group, data=samples(d))
print(design_full)

BPPARAM = MulticoreParam(ncores)

d <- dmPrecision(d, design=design_full, BPPARAM=BPPARAM)
print("estimated precision")

d <- dmFit(d, design=design_full, verbose=1, BPPARAM=BPPARAM)
print("fit model")
print(head(coefficients(d), level="feature"))

design_null <- model.matrix(~ 1, data=samples(d))

#d <- dmTest(d, coef="group1")
d <- dmTest(d, design=design_null, BPPARAM=BPPARAM)
print("finished test")

res <- results(d, level="feature")
print(head(res))

write.table(res, args$outfile, sep="\t", row.names=F, quote=F)

