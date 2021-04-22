#!/usr/local/bin/Rscript
suppressWarnings(suppressPackageStartupMessages(library(argparse)))
suppressWarnings(suppressPackageStartupMessages(library(tximport)))
suppressWarnings(suppressPackageStartupMessages(library(BANDITS)))

parser <- ArgumentParser()

parser$add_argument('--gtf', required=FALSE)
parser$add_argument('--cdna', required=FALSE)
parser$add_argument('--tx2gene', required=FALSE)
# parser$add_argument('--tool', default="kallisto", choices=c("kallisto", "salmon")) # kallisto or salmon
parser$add_argument('--pdata')
parser$add_argument('--basedir')
parser$add_argument('--ncores', default=4)
parser$add_argument('--seed', default=1984)
parser$add_argument('--outfile')

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

ncores <- as.numeric(args$ncores)

set.seed(as.numeric(args$seed))

print(args)

if(!is.null(args$tx2gene))
{
    load(args$tx2gene)
    gene_tr_id <- data.frame(gene_id=tx2gene$gene_id, transcript_id=rownames(tx2gene))

} else if(!is.null(args$gtf))
{
    suppressWarnings(suppressPackageStartupMessages(library(GenomicFeatures)))
    gtf_file <- args$gtf
    message(sprintf("reading GTF: %s...", gtf_file))
    tx <- makeTxDbFromGFF(gtf_file)
    ss <- unlist(transcriptsBy(tx, by="gene"))
    gene_tr_id <- data.frame(gene_id = names(ss), transcript_id=ss$tx_name)

} else if(!is.null(args$cdna))
{
    suppressWarnings(suppressPackageStartupMessages(library(Biostrings)))
    fasta <- readDNAStringSet(args$cdna)
    ss <- strsplit(names(fasta), " ")
    gene_tr_id <- data.frame(gene_id = gsub("gene:", "", sapply(ss, .subset, 4)),
                             transcript_id = sapply(ss, .subset, 1))

} else
{
    message("specify either --tx2gene, --cdna or --gtf to generate a transcript to gene mapping (none given)")
    q()
}

gene_tr_id <- gene_tr_id[rowSums(is.na(gene_tr_id)) == 0, ]
gene_tr_id <- unique(gene_tr_id)

print(sprintf("%d transcripts mapped to %d genes", nrow(gene_tr_id), length(unique(gene_tr_id$gene_id))))

p.data <- read.csv(args$pdata, sep="\t")
p.data$condition <- factor(p.data$condition)

print(p.data)

samples <- p.data$sample

##################################################################################
#if(args$tool == "kallisto")
#{
    message(sprintf("reading salmon input from %s for:", args$basedir))
    print(samples)

    equiv_classes <- file.path(args$basedir, samples, "aux_info", "eq_classes.txt")
    names(equiv_classes) <- samples

    if(!all(file.exists(equiv_classes)))
    {
        message("not all salmon ouput files exist!")
        print(equiv_classes)
        print(file.exists(equiv_classes))
        q()
    }
    
    quant_files = file.path(args$basedir, samples, "quant.sf")
    file.exists(quant_files)

    txi = tximport(files = quant_files, type = "salmon", txOut = TRUE)

#} else
#{
#    message("only kallisto is tool is implemented yet!")
#    q()
#}
##################################################################################

counts = txi$counts
head(counts)
eff_len = eff_len_compute(x_eff_len = txi$length)


# transcript filtering
transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
                                         transcript_counts = counts, 
                                         min_transcript_proportion = 0.01,
                                         min_transcript_counts = 10, 
                                         min_gene_counts = 20)

##################################################################################
# read equivalence class data
#if(args$tool == "kallisto")
#{
    print("reading input data from kallisto...")
    system.time(input_data <- create_data(salmon_or_kallisto="salmon",
                              gene_to_transcript=gene_tr_id,
                              salmon_path_to_eq_classes=equiv_classes,
                              eff_len=eff_len,n_cores=ncores,
                              transcripts_to_keep=transcripts_to_keep))

    print("finished reading input!")
    input_data <- filter_genes(input_data, min_counts_per_gene=20)

#} else
#{
#    message("only kallisto is tool is implemented yet!")
#    q()
#}
##################################################################################

#compute prior -> left out because of complications as reported by author
#message("computing informative prior...")
#system.time(precision <- prior_precision(gene_to_transcript=gene_tr_id,
#                             transcript_counts=counts, n_cores=ncores,
#                             transcripts_to_keep=transcripts_to_keep))
#message("finished computing prior!")

message("doing differential test...")
print(input_data)
#print(precision$prior)
print(p.data)

# save(input_data, precision, p.data, gene_tr_id, file="./input_data.RData")
# save(input_data, p.data, gene_tr_id, file="./input_data.RData")

system.time(results <- test_DTU(BANDITS_data=input_data,
                    #precision=precision$prior,
                    samples_design=p.data,
                    group_col_name="condition",
                    n_cores=ncores,
                    gene_to_transcript=gene_tr_id))
message("finished the computation")

write.table(results@Gene_results, paste0(args$outfile, ".gene.results"), sep="\t", quote=F)
write.table(results@Transcript_results, paste0(args$outfile, ".transcript.results"), sep="\t", quote=F)

