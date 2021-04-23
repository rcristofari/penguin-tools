import argparse

parser = argparse.ArgumentParser(description='Extract CDS sequences from a reference genome based on a GFF annotation file. All CDS are phased and concatenated in the order of the corresponding mRNA (not tested on anything else than the Emperor penguin reference genome)')
parser.add_argument('--ref', help='Reference genome sequence in fasta format (gzipped or not)')
parser.add_argument('--gff', help='Annotation file in GFF format (gzipped or not)')
parser.add_argument('--out', help='Basename of the output fasta file', default="extracted_cds")
parser.add_argument('--type', help='Type of output sequence (DNA or AA)', default="DNA")
parser.add_argument('--verbose', action='store_true', help="Display sequences on screen as it goes")
parser.add_argument('--legacy', action='store_true', help="For use with biopython versions < 1.78")

args = parser.parse_args()

# Imports for running the script:
import pandas as pd
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
if args.legacy:
    from Bio.Alphabet import generic_dna

# Load in the genome and parse as a dictionary:
print("Loading the reference genome...")
scaf, seq = [], []
with (gzip.open if args.ref.endswith(".gz") else open)(args.ref, "rt") as ifile:
    for record in SeqIO.parse(ifile, "fasta"):
        scaf.append(record.id)
        seq.append(str(record.seq))
genome = dict(zip(scaf, seq))

print("Loading the annotations...")
# Load in the annotation file as a dataframe:
gff = pd.read_csv(args.gff, delimiter="\t", header=None, names=["scaffold", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"])

# A function to extract a correctly phased and oriented CDS:
def extractCDS(feat, genome, verbose=True, phaseshift=0):
    # feat is a line from a gff dataframe
    # phaseshift shifts the phase set in the gff to a different value (for testing mainly)
    try:
        cds = genome[feat["scaffold"]][feat["start"] - 1:feat["end"]]
    except KeyError:
        return(None)

    if args.legacy:
        record = SeqRecord(Seq(cds, generic_dna))
    else:
        record = SeqRecord(Seq(cds))

    if feat["strand"] == "-":
        seq = str(record.reverse_complement().seq)
    else:
        seq = str(record.seq)

    cds = str(seq)
    phase = int(feat["phase"]) + phaseshift
    phased_start = cds[phase:]
    k = len(phased_start) % 3
    if k > 0:
        phased = phased_start[:-k]
    else:
        phased = phased_start

    if verbose:
        if feat["strand"] == "+":
            start, end = str(feat["start"]), str(feat["end"])
        else:
            start, end = str(feat["end"]), str(feat["start"])
        print(feat["attributes"] + " | strand: " + feat["strand"] + " | from " + start + " to " + end)
        codons = []
        for i in range(len(phased)):
            if i % 3 == 0:
                codon = phased[i:i + 3]
                if codon == "ATG":
                    codon = codon.lower()
                elif codon in ("TAG", "TGA", "TAA"):
                    print("***STOP codon found***")
                    codon = codon.lower()
                codons.append(codon)
        codons_out = " ".join(codons)
        if args.legacy:
            AA = str(SeqRecord(Seq(phased, generic_dna)).translate().seq)
        else:
            AA = str(SeqRecord(Seq(phased)).translate().seq)

        AA_out = "   ".join([x for x in AA])

        j = 0
        k = 100
        while k < len(codons_out):
            print(codons_out[j:k])
            print(AA_out[j:k])
            j += 100
            k += 100
        print(codons_out[j:])
        print(AA_out[j:])
        print("")
    return (phased)

# A function to apply the previous function to all CDS relating to one parent mRNA, and concatenate the output
def extractGene(parent, gff, genome, verbose=True, phaseshift=0, seqType="DNA", k=0, t=0):

    cds = gff.loc[gff["feature"] == "CDS"]
    cds = cds[cds['attributes'].str.contains(parent)]

    try:
        genome[cds["scaffold"].iloc[0]]
    except KeyError:
        print("Skipping gene " + parent + " (not in reference)")
        return(None)

    # Identify strand direction:
    strand = cds["strand"].iloc[0]

    print("Parsing gene " + parent + "\t| Strand: " + strand + " | CDS: " + str(cds.shape[0]) + "\t| " + str(k) + "/" + str(t))

    gene = ""
    if strand == "+":
        for i, x in cds.iterrows():
            gene = gene + extractCDS(x, genome, verbose=verbose)

    elif strand == "-":
        cds = cds.reindex(index=cds.index[::-1])
        for i, x in cds.iterrows():
            gene = gene + extractCDS(x, genome, verbose=verbose, phaseshift=phaseshift)
    else:
        raise ValueError("Unhandled strand information: " + strand)

    if verbose:
        print(">" + parent)
        print(gene)
        print("Protein sequence:")
        if args.legacy:
            aa = str(SeqRecord(Seq(gene, generic_dna)).translate().seq)
        else:
            aa = str(SeqRecord(Seq(gene)).translate().seq)
        print(aa)

    if seqType == "AA":
        return(aa)
    else:
        return (gene)

# A function to get all parent names from the GFF file:
def getAllParents(gff):
    attr = gff["attributes"].loc[gff["feature"]=="mRNA"]
    y = [x.split("=")[1] for x in attr]
    z = [x.split(";")[0] for x in y]
    return(z)

# Articulate everything:

parentList = getAllParents(gff)

with open(args.out + ".fa", "w") as ofile:
    for k, parent in enumerate(parentList):
        seq = extractGene(parent, gff, genome, seqType=args.type, verbose=args.verbose, k=k, t=len(parentList))
        if seq:
            ofile.write(">" + parent + "\n")
            ofile.write(seq + "\n")
