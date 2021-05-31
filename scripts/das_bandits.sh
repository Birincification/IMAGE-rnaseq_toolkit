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
LONGOPTS=pdata:,index:,out:,nthread:,log:,salmon,kallisto

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

sampledir="$out/SALMON/READS"
sampledir2="$out/SALMON/STAR"
sampledir3="$out/KALLISTO"

outfile="$out/diff_splicing_outs/BANDITS"
rindex="$index/R/tx2gene.RData"

for file in `find $sampledir -name "*eq_classes.txt.gz"`; do gunzip $file; done
for file in `find $sampledir2 -name "*eq_classes.txt.gz"`; do gunzip $file; done  

if [[ "$salmon" = "y" ]]; then
	watch pidstat -dru -hlH '>>' $log/bandits_salmon-$(date +%s).pidstat & wid=$!
	( [ -f "$outfile.salmon_reads.gene.results" ] && echo "$'\n'[INFO] [BANDITS] $outfile.salmon_reads already exists; skipping.." ) || \
		(/home/scripts/bandits/das_bandits.R --tx2gene $rindex --pdata $pdata --basedir $sampledir --outfile $outfile.salmon_reads --ncores $nthread --tool salmon)
	
	kill -15 $wid


	watch pidstat -dru -hlH '>>' $log/bandits_salmon_star-$(date +%s).pidstat & wid=$!

	( [ -f "$outfile.salmon_star.gene.results" ] && echo "$'\n'[INFO] [BANDITS] $outfile.salmon_star already exists; skipping.." ) || \
		(/home/scripts/bandits/das_bandits.R --tx2gene $rindex --pdata $pdata --basedir $sampledir2 --outfile $outfile.salmon_star --ncores $nthread --tool salmon)

	kill -15 $wid
fi

## for this need to update the rscript with additional params..
if [[ "$kallisto" = "y" ]]; then
	watch pidstat -dru -hlH '>>' $log/bandits_kallisto-$(date +%s).pidstat & wid=$!

	( [ -f "$outfile.kallisto.gene.results" ] && echo "$'\n'[INFO] [BANDITS] $outfile.kallisto already exists; skipping.." ) || \
		(/home/scripts/bandits/das_bandits.R --tx2gene $rindex --pdata $pdata --basedir $sampledir3 --outfile $outfile.kallisto --ncores $nthread --tool kallisto)

	kill -15 $wid
fi
## TODO implement flag for switch
#podman run -v $index:$index -v $out:$out -v $pdata:$pdata --rm -it hadziahmetovic/rnaseq-toolkit /home/scripts/das_bandits.R --tx2gene $rindex --pdata $pdata --basedir $sampledir2 --outfile $outfile.kallisto --ncores $nthread