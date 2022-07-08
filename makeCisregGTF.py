import argparse

parser = argparse.ArgumentParser(description="Create a GTF file containing cisregulatory element coordinates")
parser.add_argument('--gtf', help='Annotation file path (gz)gff')
parser.add_argument('--ref', help='Reference genome file path (gz)gff')
parser.add_argument('--kb', help='Size of the upstream cis-regulatory region, in kb [2]')
args = parser.parse_args()

import struct, os, gzip, time, datetime

start_time = time.time()

print("---------------------------------------------------------------------------------")

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
if args.gtf:
    gtf_file = args.gtf
    if not os.path.isfile(gtf_file):
        print(f"Genome file {gtf_file} could not be found")
        quit()
else:
    print("You must specify a GTF file - exiting.")
    quit()

#----------------------------------------------------------------------------------------------------------------------#
if args.kb:
    try:
        kb = int(args.kb)
        print(f"Defining cisregulatory elements as the {kb} kb area upstream of CDS")
    except ValueError:
        print(f"Invalid cisreg element size {kb}")
        quit()

print("---------------------------------------------------------------------------------")

scaf, length = [], []
scafname = None
with (gzip.open if genome_file.endswith("gz") else open)(genome_file, 'rt') as genomefile:
    for line in genomefile:
        if line.startswith(">"):
            if scafname:
                scaf.append(scafname)
                length.append(scaflength)
            scafname = line.strip(">").split(" ")[0]
            scaflength = 0
        else:
            scaflength += len(line.strip("\n"))
    scaf.append(scafname)
    length.append(scaflength)
len_dict = dict(zip(scaf, length))

if gtf_file.endswith("gz"):
    out_path = gtf_file.strip(".gtf.gz") + ".cisreg.gtf"
else:
    out_path = gtf_file.strip(".gtf") + ".cisreg.gtf"

with (gzip.open if gtf_file.endswith("gz") else open)(gtf_file, 'rt') as ifile, open(out_path, 'w') as ofile:
    for line in ifile:
        if not line.startswith("#"):
            row = line.split("\t")
            if row[2] == "start_codon":
                if row[6] == "+":
                    stop = int(row[3])
                    start = max(stop - (kb * 1000), 1)
                else:
                    start = int(row[4])
                    stop = min(start + (kb * 1000), len_dict[row[0]])
                row_out = row[:2] + ["cisreg"] + [str(start), str(stop)] + row[5:]
                ofile.write("\t".join(row_out))


#VULM01000018.1	Genbank	gene	52	1464	.	+	.	gene_id "FQA23_0006379"; transcript_id ""; gbkey "Gene"; gene_biotype "protein_coding"; locus_tag "FQA23_0006379"; partial "true";