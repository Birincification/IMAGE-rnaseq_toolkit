#!/usr/local/bin/Rscript
#args <- commandArgs(trailingOnly=T)
#if(length(args) == 0)
#{
#    print("usage: sample.file, count.dir, comp.table, gtf.file, out.dir")
#    q()
#}

#sample.file <- args[1]
#count.dir <- args[2]
#comp.table <- args[3]
#gtf.file <- args[4]
#out.dir <- args[5]
#method <- args[6]

#sample.table <- read.csv(sample.file, sep="\t", row.names="sample")
#count.files <- list.files(count.dir, full.names=T)
#comp.table <- read.csv(comp.table, sep="\t", header=F)

library(argparse)
suppressPackageStartupMessages(library(DEXSeq))


parser <- ArgumentParser()

parser$add_argument('--pdata', required=TRUE)
parser$add_argument('--countdir', required=TRUE)
parser$add_argument('--condpairs', required=TRUE)
parser$add_argument('--gtf', required=TRUE)

parser$add_argument('--outdir', required=TRUE)
parser$add_argument('--method', required=TRUE)

parser$add_argument('--ncores', default=4)
parser$add_argument('--seed', default=1984)


args <- parser$parse_args(commandArgs(trailingOnly=T))
SEED <- as.numeric(args$seed)
set.seed(SEED)
ncores <- as.numeric(args$ncores)


sample.table <- read.csv(args$pdata, sep="\t", row.names="sample")
count.files <- list.files(args$countdir, full.names=T)
comp.table <- read.csv(args$condpairs, sep="\t", header=F)
gtf.file <- args$gtf
out.dir <- args$outdir
method <- args$method

print(comp.table)
print(sample.table)
print(count.files)


tmp <- apply(comp.table, 1, function(row)
      {

        if(nrow(comp.table) > 1)
        {
            sample.table.local <- as.data.frame(sample.table[sample.table$condition==row[1] | sample.table$condition==row[2], ])
            colnames(sample.table.local) <- colnames(sample.table)
            rownames(sample.table.local) <- rownames(sample.table)[sample.table$condition==row[1] | sample.table$condition==row[2]]

            #error for most setups
            count.files.local <- count.files[c(grep(row[1], count.files),  grep(row[2], count.files))]
            #count.files.local <- count.files
        } else
        {
            sample.table.local <- sample.table
            count.files.local <- count.files
            sample.table.local$condition <- factor(sample.table.local$condition)
        }

        if(nrow(comp.table) > 1)
        {
            out.file <-   file.path(out.dir, paste0("DEXSeq.", row[1], ".vs.", row[2], ".out"))
        } else
        {
            out.file <-   file.path(out.dir, paste("DEXSeq", method, "out", sep="."))
        }


        #print(count.files)
        #print(row[1])
        #print(row[2])
        #print(sample.table.local)

        print("-----")
        print(count.files.local)
        print(sample.table.local)
        print(gtf.file)
        print(out.file)

        DDS = DEXSeqDataSetFromHTSeq(
                                     count.files.local,
                                     sampleData=sample.table.local,
                                     design= ~ sample + exon + condition:exon,
                                     flattenedfile=gtf.file )
        
        print(out.file)
        print("estimating size factors...")
        DDS <- estimateSizeFactors(DDS)
        print("estimating dispersions")

        #BiocParallel::register(BiocParallel::MulticoreParam(workers=4))
        #DDS <- estimateDispersions(DDS, BPPARAM=BiocParallel::MulticoreParam(workers=4))

		BPPARAM = MulticoreParam(ncores)

        DDS <- estimateDispersions(DDS, BPPARAM=BPPARAM)
        print("assessing DEU")
        DDS <- testForDEU(DDS, BPPARAM=BPPARAM)
        print("estimating fold changes")
        DDS = estimateExonFoldChanges( DDS, fitExpToVar="condition", BPPARAM=BPPARAM)
        print(DDS)
        dxr1 = DEXSeqResults( DDS )
        print(row)

		qvals <- perGeneQValue(dxr1)
		dxr1@listData[["qvalue"]] <- unname(ifelse(is.na(qvals[dxr1[["groupID"]]]), "1", qvals[dxr1[["groupID"]]]))

		save(dxr1, file = file.path(out.dir, "DEXSeq.object"))

        #dxr1$transcripts <- gsub("\n", "", dxr1$transcripts)
        write.table(dxr1,out.file, sep="\t", quote=F)
      })


#DDS <- DEXSeqDataSetFromSE( SE, design= ~ condition)
#DDS <- DEXSeqDataSetFromSE(SE, design= ~ sample + exon + condition:exon)

## https://bioconductor.org/packages/devel/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html