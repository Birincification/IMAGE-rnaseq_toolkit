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
LONGOPTS=pdata:,gtf:,out:,nthread:,log:

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

pdata=- out=- gtf=- nthread=4
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

gtf=$3
GTFNAME=`basename $gtf`

mkdir -p $out

#extract events from gtf
$SUPPA generateEvents -i $gtf -o $out/SUPPA2/$GTFNAME -f ioi

for d in "READS" "STAR"; do
	root=$out/SALMON/$d
	for cond in `sed '1d' $pdata | cut -f2 | sort -u`; do
		files=`grep $cond $pdata | cut -f1`
		
		## salmon is -k1 -f4; for kallisto -k1 -f5
		/home/software/SUPPA//multipleFieldSelection.py \
			-i `for file in $files; do echo $root/$file/quant.sf` -k 1 -f 4 -o $out/SUPPA2/salmon_$d.$cond.tmp.counts
		
		$SUPPA psiPerIsoform -g $gtf -e $out/SUPPA2/salmon_$(basename $root).$cond.tmp.counts -o $out/SUPPA2/$d.$cond &> /dev/null
	done

#python3.4 suppa.py joinFiles -f tpm -i sample1.tpm sample2.tpm sample3.tpm -o all_samples_tpms

##need to mute these 2 calls
#$SUPPA psiPerIsoform -g $gtf -e $out/COUNTS/tpm.counts.0 -o $out/SUPPA2/c1 &> /dev/null
#$SUPPA psiPerIsoform -g $gtf -e $out/COUNTS/tpm.counts.1 -o $out/SUPPA2/c2 &> /dev/null

	$SUPPA diffSplice -m empirical --input $out/SUPPA2/$GTFNAME.ioi --psi $out/SUPPA2/$d*_isoform.psi \
    	-e $out/SUPPA2/salmon_$d* -pa -gc -c -o $out/SUPPA2/SUPPA_salmon_$d.out

#cat $outDir/SUPPA2/SUPPA.out.dpsi.temp.0 | awk '{split($1,a,";"); print a[1]"\t"$2"\t"$3}' > $outDir/SUPPA2/test_suppa.out
	cat $out/SUPPA2/SUPPA_salmon_$d.out.dpsi | awk '{split($1,a,";"); print a[1]"\t"$2"\t"$3}' > $out/diff_splicing_outs/SUPPA_salmon_$d.out
done
#mv $OUTDIR/SUPPA.out.dpsi $5
#rm -r $OUTDIR
