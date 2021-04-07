import argparse

parser = argparse.ArgumentParser(description='Updates a (fasta) reference genome with the variants contained in a VCF file. The update can insert the IUPAC ambiguity code for polymorphic sites, or the reference / alternative / majority / minority allele. An ancestral sequence can be used to set the bases considered ambiguous in the reference sequence.')
parser.add_argument('--ref', help='Reference genome sequence in fasta format (gzipped)')
parser.add_argument('--anc', help='Ancestral genome reconstruction in fasta format (gzipped)')
parser.add_argument('--gzvcf', help='VCF file with the variants to integrate (gzipped)')
parser.add_argument('--out', help='Basename of the output fasta file', default="updated_genome")
parser.add_argument('--type', help='Which variant to use for the update (IUPAC, REF, ALT, MAJ, MIN)', default="IUPAC")
parser.add_argument('--mac', help='Minimum allele count for considering in the update', default=0)
parser.add_argument('--chrom', help='Restrict processing to this chromosome / scaffold', default=None)
parser.add_argument('--verbose', action='store_true', help="Display all runtime / debugging details")
args = parser.parse_args()

# Imports for running the script:
import gzip, re
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO

# Set some counters for final statistics:
polariseCounter, NCounter, sitesCounter, updateCounter = 0, 0, 0, 0

# Load the reference genome:
print("Loading the reference genome...")
scaf, seq = [], []
with gzip.open(args.ref, "rt") as ifile:
    for record in SeqIO.parse(ifile, "fasta"):
        scaf.append(record.id)
        seq.append(str(record.seq))
genome = dict(zip(scaf, seq))

if args.anc:
    print("Loading the ancestral genome...")
    scaf, seq = [], []
    with gzip.open(args.anc, "rt") as ifile:
        for record in SeqIO.parse(ifile, "fasta"):
            scaf.append(record.id)
            seq.append(str(record.seq))
    ancestral = dict(zip(scaf, seq))

    print("Polarising ambiguities in the reference genome...")
    reverse_iupac = {'R':['A', 'G'], 'Y':['C', 'T'], 'S':['G', 'C'],
                     'W':['A', 'T'], 'K':['G', 'T'], 'M':['A', 'C'],
                     'B':['C', 'G', 'T'], 'D':['A', 'G', 'T'],
                     'H':['A', 'C', 'T'], 'V':['A', 'C', 'G']}

    for scaf in genome:
        if args.chrom is None or scaf == args.chrom:
            refSeq = [x for x in genome[scaf]]
            for i, base in enumerate(refSeq):
                if base.upper() not in ('A', 'C', 'G', 'T', 'N'):
                    polariseCounter += 1
                    try:
                        anc_base = ancestral[scaf][i]
                        # Check that the ancestral base is compatible with the ambiguity (no double mutation on the same base)
                        if anc_base in reverse_iupac[refSeq[i].upper()]:
                            refSeq[i] = anc_base
                        else:
                            refSeq[i] = 'N'
                            NCounter += 1
                    except KeyError:
                        refSeq[i] = 'N'
                        NCounter += 1
            genome[scaf] = "".join(refSeq)
    print("Polarised the reference at " + str(polariseCounter - NCounter) + " positions")
    print("Masked " + str(NCounter) + " positions")

# Make an array for IUPAC conversion:
bases = ['A', 'C', 'G', 'T']
matrix = np.array([['A', 'M', 'R', 'W'], ['M', 'C', 'S', 'Y'], ['R', 'S', 'G', 'K'], ['W', 'Y', 'K', 'T']])
iupac = pd.DataFrame(matrix, columns=bases, index=bases)

# Start parsing the positions:
with gzip.open(args.gzvcf, 'rt') as ifile:
    for line in ifile:
        # Skip the comments
        if line[0] == "#":
            pass
        # Skip the unwanted scaffolds if the --chrom flag is set
        elif args.chrom is not None and line.split("\t")[0] != args.chrom:
            pass
        # Parse the rest
        else:
            row = line.strip("\n").split("\t")
            chrom, pos = row[0], int(row[1])
            ref, alt = row[3], row[4]
            genotypes = [x.split(':')[0] for x in row[9:]]
            alleles_pairs = [re.split("/|\|", x) for x in genotypes]
            alleles = [item for sublist in alleles_pairs for item in sublist]
            counts = Counter(alleles)
            nRef = counts['0']
            nAlt = counts['1']
            if nRef >= nAlt:
                majA, minA = ref, alt
            else:
                majA, minA = alt, ref

            if min([nRef, nAlt]) <= int(args.mac):
                ref, alt, minA = majA, majA, majA

            if args.verbose:
                print("------------------------------")
                print(chrom + ":" + str(pos) + "\t#" + str(sitesCounter + 1))
                print("REF|ALT: " + row[3] + "|" + row[4] + "\tCounts: (" + str(nRef) +"|" + str(nAlt) + ")")
                print("REF|ALT: " + ref + "|" + alt + "\tMAJ|MIN: " + majA + "|" + minA)

            #(IUPAC, REF, ALT, MAJ, MIN)
            if args.type == 'IUPAC':
                base_update = iupac[ref][alt]
            elif args.type == "REF":
                base_update = ref
            elif args.type == "ALT":
                base_update = alt
            elif args.type == "MAJ":
                base_update = majA
            elif args.type == "MIN":
                base_update = minA
            else:
                raise ValueError("Invalid base type value: " + args.type)

            if args.verbose:
                print("Old base: " + genome[chrom][pos-1] + "\tNew base: " + base_update)

            ## Insert the base at the right position in the scaffold:
            if genome[chrom][pos-1] == base_update:
                sitesCounter += 1
            else:
                genome[chrom] = genome[chrom][:pos-1] + base_update + genome[chrom][pos:]
                sitesCounter += 1
                updateCounter += 1


with open(args.out + ".fa", "w") as ofile:
    if args.chrom:
        ofile.write(">" + args.chrom + "\n")
        ofile.write(genome[args.chrom] + "\n")
    else:
        for x in genome:
            ofile.write(">" + x + "\n")
            ofile.write(genome[x] + "\n")

print("-------------------------------------------------")
print("Polarised the reference at " + str(polariseCounter) + " positions")
print("Processed " + str(sitesCounter) + " variants")
print("Updated the reference for " + str(updateCounter) + " variants")
print("-------------------------------------------------")