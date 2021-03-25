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

baseout=$out/HISAT/dta
hisat='/home/software/hisat2/hisat2'
hindex="$index/hisat2/INDEX"

dir=$(basename $out)

mkdir -p $baseout

echo "[INFO] [HISAT2] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Started processing $dir"$'\n'

for sample in `sed '1d' $pdata | cut -f1`; do
	samplein=$samples/$sample
	sampleout=$baseout/$sample.bam
	[ -f "$sampleout" ] && echo "[INFO] [HISAT] $sampleout already exists; skipping.."$'\n' && continue
	##paired
	watch pidstat -dru -hl '>>' $log/hisat2_${dir}_$sample-$(date +%s).pidstat & wid=$!

	[ -f "${samplein}_1.fastq.gz" ] &&\
		$hisat -q --dta -p $nthread -x $hindex -1 ${samplein}_1.fastq.gz -2 ${samplein}_2.fastq.gz | \
			samtools sort -@ $nthread -O BAM | tee $sampleout | samtools index -@ $nthread - $sampleout.bai
	##unpaired
	[ -f "$samplein.fastq.gz" ] &&\
		$hisat -q --dta -p $nthread -x $hindex -U $samplein.fastq.gz | \
			samtools sort -@ $nthread -O BAM | tee $sampleout | samtools index -@ $nthread - $sampleout.bai

	kill -15 $wid
done

echo "[INFO] [HISAT2] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Finished processing $dir"$'\n'
