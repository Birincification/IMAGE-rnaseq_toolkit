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
LONGOPTS=index:,pdata:,samples:,out:,nthread:,log:

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

# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
		--index)
        	index="$2"
            shift 2
            ;;
		--pdata)
            pdata="$2"
            shift 2
            ;;
        --samples)
            samples="$2"
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
        --)
            shift
            break
            ;;
        *)
            shift
            ;;
    esac
done

baseout=$out/STAR/quant
star='/home/software/STAR/bin/Linux_x86_64_static/STAR'
sindex="$index/STAR"

dir=$(basename $out)

mkdir -p $baseout

echo "[INFO] [STAR] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Started processing $dir"$'\n'

for sample in `sed '1d' $pdata | cut -f1`; do
	samplein=$samples/$sample
	sampleout=$baseout/$sample
	[ -f "$sampleout.bam" ] && echo "[INFO] [STAR] $sampleout already exists; skipping.."$'\n' && continue
	##paired
	watch pidstat -dru -hl '>>' $log/star_${dir}_$sample-$(date +%s).pidstat & wid=$!

	[ -f "${samplein}_1.fastq.gz" ] &&\
		$star --genomeDir $sindex --readFilesIn ${samplein}_1.fastq.gz ${samplein}_2.fastq.gz --runThreadN $nthread\
			--outFileNamePrefix $sampleout --readFilesCommand zcat \
			--outSAMstrandField intronMotif --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate

		samtools index -@ $nthread $sampleout.bam $sampleout.bam.bai
	##unpaired
	[ -f "$samplein.fastq.gz" ] &&\
		$star --genomeDir $sindex --readFilesIn ${samplein}.fastq.gz --runThreadN $nthread\
			--outFileNamePrefix $sampleout --readFilesCommand zcat \
			--outSAMstrandField intronMotif --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate

		samtools index -@ $nthread $sampleout.bam $sampleout.bam.bai

	kill -15 $wid
done

echo "[INFO] [STAR] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Finished processing $dir"$'\n'
