import argparse

parser = argparse.ArgumentParser(description="Create minimal per-gene GTF files containing cisregulatory element coordinates from a folder of per-gene GTFs")
parser.add_argument('--gtfd', help='Annotation file directory')
parser.add_argument('--ref', help='Reference genome file path')
args = parser.parse_args()

import struct, os, gzip, time, datetime
import pandas as pd
pd.set_option('display.max_columns', 10)
start_time = time.time()


#----------------------------------------------------------------------------------------------------------------------#
if args.ref:
    genome_file = args.ref
    if not os.path.isfile(genome_file):
        print(f"Genome file {genome_file} could not be found")
        quit()
else:
    print("You must specify a reference genome file - exiting.")
    quit()

#----------------------------------------------------------------------------------------------------------------------#
if args.gtfd:
    gtf_dir = args.gtfd
    if not os.path.isdir(gtf_dir):
        print(f"Annotation directory {gtf_file} could not be found")
        quit()
else:
    print("You must specify a GTF annotation directory file - exiting.")
    quit()

print("---------------------------------------------------------------------------------")

# read the fai file instead
scaf, length = [], []
with open(genome_file + ".fai", 'r') as indexfile:
    for line in indexfile:
        row = line.split("\t")
        scaf.append(row[0])
        length.append(int(row[1]))
len_dict = dict(zip(scaf, length))

# ADD "GENE REGION" FROM START OF CISRG TO END OF 3'UTR FOR QUICK MINING

file_list = os.listdir(gtf_dir)

for i, file in enumerate(file_list):
print(f"Processed {i} / {len(file_list)} transcripts\r", end="")
    data = pd.read_csv(gtf_dir + "/" + file, sep="\t", header=None)
    data.columns = ["chrom", "origin", "feature", "start", "end", "qual", "strand", "frame", "annotation"]

    if data[data["feature"]=="exon"].shape[0]:
        if data["strand"].iloc[0] == "+":
            ds = data[data["feature"] == "5'-UTR"]
            if ds.shape[0]:
                TSS = min(ds["start"])
                cisreg_start = max(1, TSS - 1500)
                cisreg_end = min(TSS + 500, len_dict[ds["chrom"].iloc[0]])
            else:
                ds = data[data["feature"] == "exon"]
                TSS = min(ds["start"])
                cisreg_start = max(1, TSS - 2000)
                cisreg_end = TSS

        else:
            ds = data[data["feature"] == "5'-UTR"]
            if ds.shape[0]:
                TSS = max(ds["end"])
                cisreg_start = max(1, TSS - 500)
                cisreg_end = min(TSS + 1500, len_dict[ds["chrom"].iloc[0]])
            else:
                ds = data[data["feature"] == "exon"]
                TSS = max(ds["end"])
                cisreg_start = TSS
                cisreg_end = min(TSS + 2000, len_dict[ds["chrom"].iloc[0]])

        cisreg = {"chrom": data["chrom"][0], "origin": "python3", "feature": "cisreg", "start": cisreg_start,
                 "end": cisreg_end, "qual": ".", "strand": data["strand"][0], "frame": ".",
                 "annotation": data["annotation"][0]}

        data = pd.concat([data, pd.DataFrame.from_records([cisreg])], ignore_index=True)
        data = data.sort_values(["start", "end"])

        data.to_csv(gtf_dir + "/" + file.strip(".gtf") + ".cisreg.gtf", sep='\t', header=False, index=False)
