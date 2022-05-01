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

diffexp_script="/home/scripts/deseq/de_rseq.R"
create_input='/home/scripts/deseq/create_exprs_fdata.py'

samples=()
for name in `sed '1d' $pdata | cut -f1 | awk '{print $1 ".bam"}'`; do
        samples+=( "$name" )
done


mkdir -p $out/diff_exp_outs
#DESeq
for method in "hisat" "star" "contextmap" "ideal"; do
	if [[ "${map[$method]}" = "y" ]]; then

		basein=$out/COUNTS/featureCounts.$method

		watch pidstat -dru -hlH '>>' $log/deseq_$method-$(date +%s).pidstat & wid=$!

		##generate separate count files
		echo "[INFO] [DiffExp] ["`date "+%Y/%m/%d-%H:%M:%S"`"] $create_input"
		$create_input $basein $method $pdata $out/COUNTS ${samples[@]}

		echo "[INFO] [DiffExp] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Preprocessing finished"$'\n'

		kill -15 $wid


		watch pidstat -dru -hlH '>>' $log/deseq_$method-$(date +%s).pidstat & wid=$!

		( [ -f "$out/diff_exp_outs/DESeq.out" ] && echo "[INFO] [DESeq] $out/diff_exp_outs/DESeq.out already exists, skipping.."$'\n' ) || \
            ( $diffexp_script $out/COUNTS/exprs_$method.txt $out/COUNTS/p_data_$method.txt $out/COUNTS/f_data_$method.txt DESeq $out/diff_exp_outs/DESeq_$method.out )
 

        ( [ -f "$out/diff_exp_outs/edgeR.out" ] && echo "[INFO] [edgeR] $out/diff_exp_outs/edgeR.out already exists, skipping.."$'\n' ) || \
            ( $diffexp_script $out/COUNTS/exprs_$method.txt $out/COUNTS/p_data_$method.txt $out/COUNTS/f_data_$method.txt edgeR $out/diff_exp_outs/edgeR_$method.out )


        ( [ -f "$out/diff_exp_outs/limma.out" ] && echo "[INFO] [limma] $out/diff_exp_outs/limma.out already exists, skipping.."$'\n' ) || \
            ( $diffexp_script $out/COUNTS/exprs_$method.txt $out/COUNTS/p_data_$method.txt $out/COUNTS/f_data_$method.txt limma $out/diff_exp_outs/limma_$method.out )

		kill -15 $wid
	fi
done

#[INFO] [edgeR] [2022/05/01-12:06:53] /home/scripts/de_rseq_new.R /home/data/out//COUNTS/exprs.txt out/COUNTS/p_data.txt /home/data/out//COUNTS/f_data.txt edgeR /home/data/out//diff_exp_outs/edgeR.out

echo "[INFO] [DiffExp] ["`date "+%Y/%m/%d-%H:%M:%S"`"] DiffExp analysis finished"$'\n'

#        ( [ -f "$out/diff_exp_outs/sleuth.out" ] && echo "[INFO] [Sleuth] $out/diff_exp_outs/sleuth.out already exists, skipping.."$'\n' ) || \
#            ( /home/scripts/run_sleuth.R --base.dir $out/KALLISTO/alignment --p.data $pdata --condition.pairs /home/data/cond.pairs \
#                --out.dir $out/diff_exp_outs --num.cores $nthread --tx2gene $index/R/tx2gene.RData )
#
#
#        ( [ -f "$out/diff_exp_outs/ballgown.out" ] && echo "[INFO] [Ballgown] $out/diff_exp_outs/ballgown.out already exists, skipping.."$'\n' ) || \
#            ( /home/scripts/run_ballgown.R --base.dir $out/STRINGTIE --p.data $pdata --condition.pairs /home/data/cond.pairs \
#                --out.dir $out/diff_exp_outs --tx2gene $index/R/tx2gene.RData )
