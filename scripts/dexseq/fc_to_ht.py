#!/usr/bin/python3
import argparse
import pandas as pd
import os

def main(args):
    #sample_table = pd.read_csv(args.sampletable, sep="\t")
    #print(sample_table)

    with open(args.sampletable, 'r') as infile:
        sample_table = [l[:-1] for l in infile.readlines()]


    if not os.path.exists(args.htdir):
        os.makedirs(args.htdir)

    #sample_files = {sample:open(os.path.join(args.htdir, sample), "w")  for sample in sample_table["sample"]}
    sample_files = {sample:open(os.path.join(args.htdir, sample), "w")  for sample in sample_table}


    counts = pd.read_csv(args.fcfile, sep="\t", comment="#", dtype={"Chr" : str})
    
    for _, row in counts.iterrows():
        
        tmp = [float(row[x]) for x in sample_files.keys()]

        if sum(tmp) == 0:
            continue


        for sample, sample_file in sample_files.items():
            #if row[sample] == 0:
            #    continue
            sample_file.write("%s\t%d\n"%(row["Geneid"], row[sample]))



    for f in sample_files.values():
        f.write("_ambiguous\t0\n_ambiguous_readpair_position\t0\n_empty\t0\n_lowaqual\t0\n_notaligned\t0\n")
        f.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--fcfile", required=True)
    parser.add_argument("--htdir", required=True)
    parser.add_argument("--sampletable", required=True)

    args = parser.parse_args()

    main(args)
