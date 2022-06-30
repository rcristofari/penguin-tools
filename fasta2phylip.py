import argparse

parser = argparse.ArgumentParser(description='Takes a tab-separated file of species codes (max 4 chars) and genome file paths, and produces scaffold-by-scaffold Phylip alignments, for use with baseml.')
parser.add_argument('--files', help='Tab-separated species names / file path list')
args = parser.parse_args()

from Bio import SeqIO
import time

species, paths = [], []
with open(args.files, 'r') as ifile:
    for line in ifile:
        spe, path = line.strip("\n").split("\t")
        species.append(spe)
        paths.append(path)

fastas = []
for path in paths:
    fastas.append(list(SeqIO.parse(path, "fasta")))

n_records = len(fastas[0])
n_species = len(species)

print(f"Processing data for {n_records} scaffolds")

for i in range(n_records):
    file_name = fastas[0][i].id + ".phy"
    with open(file_name, 'w') as ofile:
        ofile.write(f"{n_species}\t{len(fastas[0][i].seq)}\n")
        for j, s in enumerate(species):
            ofile.write(s.ljust(5, " "))
            ofile.write(str(fastas[j][i].seq))
            ofile.write("\n")
    if i % 100 == 0:
        print(f"{time.asctime()} | Processed {i} scaffolds", end="\r")