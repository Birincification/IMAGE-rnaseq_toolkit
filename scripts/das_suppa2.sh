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
LONGOPTS=pdata:,gtf:,out:,nthread:,log:,salmon,salmonstar,kallisto

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

pdata=- out=- gtf=- nthread=4 salmon=n salmonstar=n kallisto=n
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        --pdata)
            pdata="$2"
            shift 2
            ;;
		--gtf)
        	gtf="$2"
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
        --kallisto)
            kallisto=y
            shift
            ;;
        --salmon)
            salmon=y
            shift
            ;;
        --salmonstar)
            salmonstar=y
            shift
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

# handle non-option arguments
if [[ $# -ne 0 ]]; then
    echo "$0: empty flag detected, is this intentional?!"
    #exit 4
fi

SUPPA="python3 /home/software/SUPPA/suppa.py"

GTFNAME=`basename $gtf`

dir=$(basename $out)

mkdir -p $out/SUPPA2

#extract events from gtf
$SUPPA generateEvents -i $gtf -o $out/SUPPA2/$GTFNAME -f ioi

if [[ "$salmon" = "y" ]]; then
	d="READS"
	root=$out/SALMON/$d

	watch pidstat -dru -hlH '>>' $log/suppa2-psiPeform_${dir}_salmon.$(date +%s).pidstat & wid=$!
	starter="$(date +%s)"

	for cond in `sed '1d' $pdata | cut -f2 | sort -u`; do
		files=`grep -P "\t$cond" $pdata | cut -f1`
		f=`for file in $files; do echo $root/$file/quant.sf; done`
		

		## salmon is -k1 -f4; for kallisto -k1 -f5
		python3 /home/software/SUPPA//multipleFieldSelection.py \
			-i $f -k 1 -f 4 -o $out/SUPPA2/salmon_$cond.$d.tmp.counts

		python3 /home/software/SUPPA//multipleFieldSelection.py \
			-i $f -k 1 -f 5 -o $out/SUPPA2/salmon_$cond.$d.read.counts
		
		$SUPPA psiPerIsoform -g $gtf -e $out/SUPPA2/salmon_$cond.$d.tmp.counts -o $out/SUPPA2/$d.$cond
	done
	echo "$(($(date +%s)-$starter))" >> $log/suppa2-psiPeform_${dir}_salmon.$(date +%s).runtime
	kill -15 $wid


	watch pidstat -dru -hlH '>>' $log/suppa2-diff_${dir}_salmon.$(date +%s).pidstat & wid=$!
	starter="$(date +%s)"

	( [ -f "$out/diff_splicing_outs/SUPPA_salmon_$d.out" ] && echo "[INFO] [SUPPA2] $out/diff_splicing_outs/SUPPA_salmon_$d.out already exists, skipping.."$'\n' ) \
			|| ($SUPPA diffSplice -m empirical --input $out/SUPPA2/$GTFNAME.ioi --psi $out/SUPPA2/$d*_isoform.psi \
    	-e $out/SUPPA2/salmon_*.$d.tmp.counts -gc -o $out/SUPPA2/SUPPA_salmon_$d.out)

	echo "$(($(date +%s)-$starter))" >> $log/suppa2-diff_${dir}_salmon.$(date +%s).runtime
	kill -15 $wid

	( [ -f "$out/diff_splicing_outs/SUPPA_salmon_$d.out" ]) || \
	(cat $out/SUPPA2/SUPPA_salmon_$d.out.dpsi | awk '{split($1,a,";"); print a[1]"\t"$2"\t"$3}' | sort -r -u > $out/diff_splicing_outs/SUPPA_salmon_$d.out)
fi


if [[ "$salmonstar" = "y" ]]; then
	d="STAR"
	root=$out/SALMON/$d

	watch pidstat -dru -hlH '>>' $log/suppa2-psiPeform_${dir}_salmon-star.$(date +%s).pidstat & wid=$!
	starter="$(date +%s)"

	for cond in `sed '1d' $pdata | cut -f2 | sort -u`; do
		files=`grep -P "\t$cond" $pdata | cut -f1`
		f=`for file in $files; do echo $root/$file/quant.sf; done`
		
		## salmon is -k1 -f4; for kallisto -k1 -f5
		python3 /home/software/SUPPA//multipleFieldSelection.py \
			-i $f -k 1 -f 4 -o $out/SUPPA2/salmon_$cond.$d.tmp.counts

		python3 /home/software/SUPPA//multipleFieldSelection.py \
			-i $f -k 1 -f 5 -o $out/SUPPA2/salmon_$cond.$d.read.counts
		
		$SUPPA psiPerIsoform -g $gtf -e $out/SUPPA2/salmon_$cond.$d.tmp.counts -o $out/SUPPA2/$d.$cond &> /dev/null
	done
	echo "$(($(date +%s)-$starter))" >> $log/suppa2-psiPeform_${dir}_salmon-star.$(date +%s).runtime
	kill -15 $wid


	watch pidstat -dru -hlH '>>' $log/suppa2-diff_${dir}_salmon-star.$(date +%s).pidstat & wid=$!
	starter="$(date +%s)"

	( [ -f "$out/diff_splicing_outs/SUPPA_salmon_$d.out" ] && echo "[INFO] [SUPPA2] $out/diff_splicing_outs/SUPPA_salmon_$d.out already exists, skipping.."$'\n' ) \
			|| ($SUPPA diffSplice -m empirical --input $out/SUPPA2/$GTFNAME.ioi --psi $out/SUPPA2/$d*_isoform.psi \
    	-e $out/SUPPA2/salmon_*.$d.tmp.counts -gc -o $out/SUPPA2/SUPPA_salmon_$d.out)

	echo "$(($(date +%s)-$starter))" >> $log/suppa2-diff_${dir}_salmon-star.$(date +%s).runtime
	kill -15 $wid

	( [ -f "$out/diff_splicing_outs/SUPPA_salmon_$d.out" ]) || \
	(cat $out/SUPPA2/SUPPA_salmon_$d.out.dpsi | awk '{split($1,a,";"); print a[1]"\t"$2"\t"$3}' | sort -r -u > $out/diff_splicing_outs/SUPPA_salmon_$d.out)
fi


if [[ "$kallisto" = "y" ]]; then
	d="quant"
	root=$out/KALLISTO/$d

	watch pidstat -dru -hlH '>>' $log/suppa2-psiPeform_${dir}_kallisto.$(date +%s).pidstat & wid=$!
	starter="$(date +%s)"

	for cond in `sed '1d' $pdata | cut -f2 | sort -u`; do
		files=`grep -P "\t$cond" $pdata | cut -f1`
		f=`for file in $files; do echo $root/$file/abundance.tsv; done`
		

		## salmon is -k1 -f4; for kallisto -k1 -f5
		python3 /home/software/SUPPA//multipleFieldSelection.py \
			-i $f -k 1 -f 5 -o $out/SUPPA2/kallisto_$cond.tmp.counts

		python3 /home/software/SUPPA//multipleFieldSelection.py \
			-i $f -k 1 -f 4 -o $out/SUPPA2/kallisto_$cond.read.counts
		
		$SUPPA psiPerIsoform -g $gtf -e $out/SUPPA2/kallisto_$cond.tmp.counts -o $out/SUPPA2/kallisto.$cond &> /dev/null
	done
	echo "$(($(date +%s)-$starter))" >> $log/suppa2-psiPeform_${dir}_kallisto.$(date +%s).runtime
	kill -15 $wid

	watch pidstat -dru -hlH '>>' $log/suppa2-diff_${dir}_kallisto.$(date +%s).pidstat & wid=$!
	starter="$(date +%s)"

	( [ -f "$out/diff_splicing_outs/SUPPA_kallisto.out" ] && echo "[INFO] [SUPPA2] $out/diff_splicing_outs/SUPPA_kallisto.out already exists, skipping.."$'\n' ) \
			|| ($SUPPA diffSplice -m empirical --input $out/SUPPA2/$GTFNAME.ioi --psi $out/SUPPA2/kallisto*_isoform.psi \
    	-e $out/SUPPA2/kallisto*.tmp.counts -gc -o $out/SUPPA2/SUPPA_kallisto.out)

	echo "$(($(date +%s)-$starter))" >> $log/suppa2-diff_${dir}_kallisto.$(date +%s).runtime
	kill -15 $wid

	( [ -f "$out/diff_splicing_outs/SUPPA_kallisto.out" ]) || \
	(cat $out/SUPPA2/SUPPA_kallisto.out.dpsi | awk '{split($1,a,";"); print a[1]"\t"$2"\t"$3}' | sort -r -u > $out/diff_splicing_outs/SUPPA_kallisto.out)
fi