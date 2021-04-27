#!/usr/bin/python3
import sys
infile = sys.argv[1]
outfile = infile+".filtered"


with open(infile, "r") as incount, open(outfile, "w") as outcount:
    current_gene=""
    i = 1

    for l in incount:
        if l.startswith("#"):
            outcount.write(l)
            continue


        #if l.startswith("Geneid"):
        #    outcount.write(l)
        #    continue

        sp = l.split("\t")


        if l.startswith("Geneid"):
            for ii, ss in enumerate(sp):
                if not "/" in ss:
                    continue

                sp[ii] = ".".join(sp[ii].split("/")[-1].split(".")[:-1])

                #sp2 = ss.split("/")
                #sp3 = sp2[-2].split("_")
                #sp[ii] = "%s.%s%s"%(sp3[0], sp3[1], sp3[2])
            outcount.write("\t".join(sp) + "\n")
            continue
            
        this_gene = sp[0]

        if not this_gene == current_gene:
            i = 1
            current_gene = this_gene

        sp[0] = "%s:%.3d"%(sp[0],i)
        outcount.write("\t".join(sp))
        i+=1


        
