#!/bin/bash -x
## -x : expands variables and prints commands as they are called

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
LONGOPTS=pdata:,index:,out:,nthread:,log:,hisat2,star,contextmap,ideal,

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

declare -A map

pdata=- out=- index=- nthread=4 map[hisat]=n map[star]=n
map[contextmap]=n map[ideal]=n
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        --pdata)
            pdata="$2"
            shift 2
            ;;
		--index)
        	index="$2"
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
	    --hisat2)
            map[hisat]=y
            shift
            ;;
        --star)
            map[star]=y
            shift
            ;;
        --contextmap)
            map[contextmap]=y
            shift
            ;;
        --ideal)
            map[ideal]=y
            shift
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

# handle non-option arguments
if [[ $# -ne 0 ]]; then
    echo "$0: empty flag detected, is this intentional?!"
    #exit 4
fi


sed '1d' $pdata | cut -f2 | uniq | tr "\\n" "\\t" | rev | cut -c 2- | rev > /home/cond.pairs
sed '1d' $pdata | cut -f1 | awk '{print $1 ".bam"}' > /home/sample.list

dexseq_script="/home/scripts/dexseq/das_dexseq.R"
dexseqFilter='/home/scripts/dexseq/modify_dexseq_featurecounts.py'
dexseqHT='/home/scripts/dexseq/fc_to_ht.py'

dir=$(basename $out)

mkdir -p $out/DEXSEQ
#DEXSeq
for method in "hisat" "star" "contextmap" "ideal"; do
	if [[ "${map[$method]}" = "y" ]]; then

		basein=$out/COUNTS/featureCounts.$method.DEXSeq

		watch pidstat -dru -hlH '>>' $log/dexseqFilter_${dir}_$method.$(date +%s).pidstat & wid=$!

		##generate counts
		echo "[INFO] [DEXSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] $dexseqFilter"

		## here add $method to DEXSeq.$method and also in fC generation
		$dexseqFilter $basein
		echo "[INFO] [DEXSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] $dexseqHT"
		$dexseqHT --fcfile $basein.filtered --htdir $out/DEXSEQ/${method}_HTcounts --sampletable /home/sample.list
		echo "[INFO] [DEXSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Preprocessing finished"$'\n'

		kill -15 $wid


		watch pidstat -dru -hlH '>>' $log/dexseq_${dir}_$method.$(date +%s).pidstat & wid=$!

		( [ -f "$out/diff_splicing_outs/DEXSeq.$method.out" ] && echo "[INFO] [DEXSeq] $out/diff_splicing_outs/DEXSeq.$method.out already exists, skipping.."$'\n' ) \
			|| ($dexseq_script --pdata $pdata --countdir $out/DEXSEQ/${method}_HTcounts --condpairs /home/cond.pairs \
				--gtf $index/dexseq/annot.noaggregate.gtf --outdir $out/diff_splicing_outs --method $method --ncores $nthread)

		kill -15 $wid
	fi
done

echo "[INFO] [DEXSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Splicing analysis finished"$'\n'