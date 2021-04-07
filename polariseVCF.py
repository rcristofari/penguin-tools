import argparse

parser = argparse.ArgumentParser(description='Polarises a VCF file according to a reference other than the one originally used for SNP calling (but sharing the same coordinates, typically an ancestral reconstruction)')
parser.add_argument('--ref', help='Reference genome sequence in fasta format (gzipped)')
parser.add_argument('--gzvcf', help='VCF file to polarise (gzipped)')
parser.add_argument('--out', help='Basename of the output fasta file', default="polarised")
args = parser.parse_args()

import sys, os, gzip
from Bio import SeqIO
import numpy as np

scaf, seq = [], []
with gzip.open(args.ref, "rt") as ifile:
    for record in SeqIO.parse(ifile, "fasta"):
        scaf.append(record.id)
        seq.append(str(record.seq))

genome = dict(zip(scaf, seq))

lineN, reverseN = 0, 0

with gzip.open(args.gzvcf, "rt") as ifile:
    with open(args.out + ".vcf", "w") as ofile:
        for line in ifile:

            # For header lines:
            if line[0] == "#":
                ofile.write(line)
       # For record lines:
            else:
                lineN += 1
                row = line.strip("\n").split("\t")
                s = row[0]
                p = int(row[1])
                ref = row[3]
                alt = row[4]
                anc = genome[s][p-1]


                if anc == alt:
                    reverseGT = True
                else:
                    reverseGT = False
                    ofile.write(line)

                # If we need to reverse everything:
                if reverseGT:
                    reverseN += 1
                    # Initialise the row:
                    line_out = row[:9]
                    line_out[3] = alt
                    line_out[4] = ref

                    # We then parse all the rest of the data:
                    fields = np.array(row[8].split(":"))
                    fieldDict = {"GT":None, "AD":None, "DP":None, "GQ":None, "PGT":None, "PID":None, "PL":None, "PS":None}

                    for g in row[9:]:

                        # Genotype:
                        if "GT" in fields:
                            gtloc = int(np.where(fields == "GT")[0])
                            gt = g.split(":")[gtloc]
                            if "/" in gt:
                                gt = gt.split("/")
                                gt.reverse()
                                gt = "/".join(gt)
                            elif "|" in gt:
                                gt = gt.split("|")
                                gt.reverse()
                                gt = "|".join(gt)
                            else:
                                raise ValueError("Unrecognised genotype format: " + gt)
                            fieldDict["GT"] = gt

                        # Allelic depth:
                        if "AD" in fields:
                            adloc = int(np.where(fields == "AD")[0])
                            ad = g.split(":")[adloc].split(",")
                            ad.reverse()
                            ad = ",".join(ad)
                            fieldDict["AD"] = ad

                        # Overall depth:
                        if "DP" in fields:
                            dploc = int(np.where(fields == "DP")[0])
                            dp = g.split(":")[dploc]
                            fieldDict["DP"] = dp

                        # Genotype quality:
                        if "GQ" in fields:
                            gqloc = int(np.where(fields == "GQ")[0])
                            gq = g.split(":")[gqloc]
                            fieldDict["GQ"] = gq

                        # Physical phasing:
                        if "PGT" in fields:
                            pgtloc = int(np.where(fields == "PGT")[0])
                            pgt = g.split(":")[pgtloc]
                            if "|" in pgt:
                                pgt = pgt.split("|")
                                pgt.reverse()
                                pgt = "|".join(pgt)
                            else:
                                pgt = "."
                            fieldDict["PGT"] = pgt

                        # Physical phasing ID:
                        if "PID" in fields:
                            pidloc = int(np.where(fields == "PID")[0])
                            pid = g.split(":")[pidloc]
                            fieldDict["PID"] = pid

                        # Phred-scaled likelihoods:
                        if "PL" in fields:
                            plloc = int(np.where(fields == "PL")[0])
                            pl = g.split(":")[plloc].split(",")
                            pl.reverse()
                            pl = ",".join(pl)
                            fieldDict["PL"] = pl

                        # Phasing set:
                        if "PS" in fields:
                            psloc = int(np.where(fields == "PS")[0])
                            ps = g.split(":")[psloc]
                            fieldDict["PS"] = ps

                        if any(x not in ("GT", "AD", "DP", "GQ", "PGT", "PID", "PL", "PS") for x in fields):
                            raise ValueError('Unhandled field name')


                        geno_list = [fieldDict[x] for x in fields]

                        line_out.append(":".join(geno_list))

                    ofile.write("\t".join(line_out) + "\n")

sys.stdout.write("Finished parsing %d positions\n" % (lineN) )
sys.stdout.write("Repolarised %d positions\n" % (reverseN) )
