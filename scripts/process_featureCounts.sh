#!/bin/bash -x

echo $@
params=("$@")

# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi

OPTIONS=
LONGOPTS=index:,gtf:,pdata:,out:,nthread:,log:,dexseq

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

dexseq=n
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
		--index)
        	index="$2"
            shift 2
            ;;
		--gtf)
        	gtf="$2"
            shift 2
            ;;
		--pdata)
            pdata="$2"
            shift 2
            ;;
        --out)
            out="$2"
            shift 2
            ;;
        --nthread)
            nthread="$2"
            shift 2
            ;;
		--log)
			log="$2"
			shift 2
			;;
		--dexseq)
			dexseq="y"
			shift 2
			;;
        --)
            shift
            break
            ;;
        *)
            shift
            ;;
    esac
done

basein=$out/HISAT/dta
baseout=$out/COUNTS/featureCounts
featureCounts='/home/software/subread/bin/featureCounts'

dir=$(basename $out)

mkdir -p $out/COUNTS

[ -f "$baseout" ] && echo "[INFO] [featureCounts] $baseout already exists; skipping.."$'\n' && exit 

hisatBams=()
for sample in `sed '1d' $pdata | cut -f1`; do
	hisatBams+=( "$sample.bam" )
done

echo "[INFO] [featureCounts] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Started processing $baseout"$'\n'
##gene level counting
cd $out/HISAT/dta
watch pidstat -dru -hl '>>' $log/featureCounts_${dir}-$(date +%s).pidstat & wid=$!

$featureCounts -T $nthread -B -C -a $gtf --primary -o $baseout ${hisatBams[@]}

kill -15 $wid

##exon_part level counting
if [[ "$dexseq" = "y" ]]; then
	echo "[INFO] [featureCounts] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Started processing $baseout.DEXSeq"$'\n'
	watch pidstat -dru -hl '>>' $log/featureCounts_exonicpart_${dir}-$(date +%s).pidstat & wid=$!
	
	$featureCounts -T $nthread -B -C -O -f -t exonic_part -a $index/dexseq/annot.noaggregate.gtf --primary -o $baseout.DEXSeq ${hisatBams[@]}

	kill -15 $wid
fi

echo "[INFO] [featureCounts] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Finished processing $out"$'\n'

#  -a <string>         Name of an annotation file. GTF/GFF format by default.
#                      See -F option for more format information. Inbuilt
#                      annotations (SAF format) is available in 'annotation'
#                      directory of the package.
#
#  -o <string>         Name of the output file including read counts. A separate
#                      file including summary statistics of counting results is
#                      also included in the output ('<string>.summary').
#
#  -s <int>            Perform strand-specific read counting. Acceptable values:
#                      0 (unstranded), 1 (stranded) and 2 (reversely stranded).
#                      0 by default.
#
#  -T <int>            Number of the threads. 1 by default.
#
#  -p                  If specified, fragments (or templates) will be counted
#                      instead of reads. This option is only applicable for
#                      paired-end reads.
#
#  -B                  Only count read pairs that have both ends aligned.
#
#  -C                  Do not count read pairs that have their two ends mapping 
#                      to different chromosomes or mapping to same chromosome 
#                      but on different strands.
#
#  --fracOverlap <float> Minimum fraction of overlapping bases in a read that is
#                      required for read assignment. Value should be within range
#                      [0,1]. 0 by default. Number of overlapping bases is
#                      counted from both reads if paired end. Both this option
#                      and '--minOverlap' option need to be satisfied for read
#                      assignment.
#
#  --primary           Count primary alignments only. Primary alignments are 
#                      identified using bit 0x100 in SAM/BAM FLAG field.
#
#  -t <string>         Specify feature type in GTF annotation. 'exon' by 
#                      default. Features used for read counting will be 
#                      extracted from annotation using the provided value.
#
#  -O                  Assign reads to all their overlapping meta-features (or 
#                      features if -f is specified).
#
#  -f                  Perform read counting at feature level (eg. counting 
#                      reads for exons rather than genes).
#
#/home/software/subread/bin/featureCounts -s 0 -T 6 -p -B -C --fracOverlap 1 -a /home/GENOMIC_DERIVED/PAN/homo_sapiens/standardchr/Homo_sapiens.GRCh37.75.gtf --primary -o output/COUNTS/featureCounts output/cond1_02.bam output/cond2_01.bam output/cond1_01.bam output/cond2_00.bam output/cond1_00.bam output/cond2_02.bam

#/home/software/subread/bin/featureCounts -O -s 0 -T 6 -p -B -C -f -t exonic_part -a /home/indices/DEXSeq/9606/standardchr/annot.noaggregate.gtf --primary -o output/COUNTS/featureCounts.DEXSeq output/cond1_02.bam output/cond2_01.bam output/cond1_01.bam output/cond2_00.bam output/cond1_00.bam output/cond2_02.bam
