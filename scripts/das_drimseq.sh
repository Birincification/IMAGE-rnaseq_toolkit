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
LONGOPTS=pdata:,index:,out:,nthread:,log:,salmon,kallisto,

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

pdata=- out=- index=- nthread=4 salmon=n kallisto=n
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
        --salmon)
            salmon=y
            shift
            ;;
        --kallisto)
            kallisto=y
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

drimseq="/home/scripts/drimseq/das_drimseq.R"

sampledir="$out/SALMON/READS"
sampledir2="$out/SALMON/STAR"
sampledir3="$out/KALLISTO/quant"

outfile="$out/diff_splicing_outs/DRIMSeq"
rindex="$index/R/tx2gene.RData"

#DRIMSeq on kallisto
if [[ "$kallisto" = "y" ]]; then
	watch pidstat -dru -hlH '>>' $log/drimseq_kallisto-$(date +%s).pidstat & wid=$!
	
	( [ -f "$outfile.kallisto.out" ] && echo "[INFO] [DRIMSeq] $outfile.kallisto.out already exists, skipping.."$'\n' ) \
		|| ($drimseq --counts $sampledir3 --pdata $pdata --outfile $outfile.kallisto.out --tx2gene $rindex --ncores $nthread)
	
	kill $wid
fi

echo "[INFO] [DRIMSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Splicing analysis finished"$'\n'
