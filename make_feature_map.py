import argparse, struct, os, gzip, time, datetime
import numpy as np

parser = argparse.ArgumentParser(description="Create a feature map for a reference genome from a GFF file and a CpG cluster map. Create a binary methylation matrix file referenced to genome features, for a given set of individuals.")
parser.add_argument('--ref', help='Reference genome path (gzipped fasta file)')
parser.add_argument('--gff', help='Annotation file path (gzipped gff file)')
parser.add_argument('--cpgi', help='CpG cluster file (txt file from gCluster)')
parser.add_argument('--vcf', help='SNP database in VCF format, to annotate C/T polymorphisms')
parser.add_argument('--cpg', help='Annotate CpG sites from the reference genome')
parser.add_argument('--map', help='Binary feature map file output by a previous run')
parser.add_argument('--matrix', help='Methylation matrix (one position per line, one sample per column, no missing values)')
parser.add_argument('--samples', help='Sample list (text file, one sample per row)')
parser.add_argument('--output', help='Output folder path')
args = parser.parse_args()

start_time = time.time()

genome_file = args.ref
if not os.path.isfile(genome_file):
    print(f"Genome file {genome_file} could not be found")
    quit()

if args.gff:
    doGFF = True
    genome_gff = args.gff
    if not os.path.isfile(genome_gff):
        print(f"GFF file {genome_gff} could not be found")
        quit()
else:
    doGFF = False

if args.cpgi:
    doCPGI = True
    genome_cpgi = args.cpgi
    if not os.path.isfile(genome_cpgi):
        print(f"CpG cluster file {genome_cpgi} could not be found")
        quit()
else:
    doCPGI = False

if args.vcf:
    doVCF = True
    dbsnp_path = args.vcf
    if not os.path.isfile(dbsnp_path):
        print(f"VCF SNP file {dbsnp_path} could not be found")
        quit()
else:
    doVCF = False

if args.cpg:
    doCpG = True
else:
    doCpG = False

if args.map:
    doMap = False
    map_path = args.map
    if not os.path.isfile(map_path):
        print(f"Map file {map_path} could not be found")
        quit()
else:
    doMap = True
    map_path = out_dir + "features.map"


if args.matrix:
    doMatrix = True
    matrix_file = args.matrix
    if not os.path.isfile(matrix_file):
        print(f"Matrix file {matrix_file} could not be found")
        quit()

    sample_file = args.samples
    if not os.path.isfile(sample_file):
        print(f"Sample list {sample_file} could not be found")
        quit()
    else:
        samples = []
        with open(sample_file, 'r') as ifile:
            for line in ifile:
                samples.append(line.strip("\n"))
        print(f"Processing data for {len(samples)} samples.")
else:
    doMatrix = False

if args.output:
    out_dir = args.output
    if not os.path.isdir(out_dir):
        print(f"Output directory {out_dir} could not be found, defaulting to current directory")
        out_dir = "./"
    if not out_dir.endswith("/"):
        out_dir += "/"
else:
    out_dir = "./"

#------------------------------------------------------------------------------#
# The exponents used for the bit codes:

codes = {'cpg':0, 'cpgi':1, 'cisreg':2, 'intron':3, 'exon':4, 'ct_snp':5}

#------------------------------------------------------------------------------#
# Figure out the total length of the reference genome:

print("# Calculating target genome size...", end='')
genome_size = 0
with gzip.open(genome_file, 'rt') as gfile:
    for line in gfile:
        if line.startswith(">"):
            pass
        else:
            genome_size += len(line.strip('\n'))
print(f"done. Genome size = {genome_size} bp")

#------------------------------------------------------------------------------#
# Create a coordinate reference file for the scaffolds in the genome:

print("# Indexing scaffolds...", end='')
coord_file = out_dir + "genome.coord"
cumulative_start, scafname, scaflength = 0, None, 0
scaf, length, abs_start = [], [], []

with (gzip.open if genome_file.endswith("gz") else open)(genome_file, 'rt') as genomefile, open(coord_file, 'w') as coordfile:
    for line in genomefile:
        if line.startswith(">"):
            if scafname:
                coordfile.write("\t".join([scafname, str(scaflength), str(cumulative_start)]) + "\n")
                scaf.append(scafname)
                length.append(scaflength)
                abs_start.append(cumulative_start)

            scafname = line.strip(">").split(" ")[0]
            cumulative_start += scaflength
            scaflength = 0

        else:
            scaflength += len(line.strip("\n"))
    coordfile.write("\t".join([scafname, str(scaflength), str(cumulative_start)]) + "\n")
    scaf.append(scafname)
    length.append(scaflength)
    abs_start.append(cumulative_start)

len_dict = dict(zip(scaf, length))
abs_dict = dict(zip(scaf, abs_start))

print("done.")

#------------------------------------------------------------------------------#
# Create the feature map file

if doMap:
    print("\n# Initializing the feature map file...")
    # Initialise an empty code file:
    with open(map_path, 'wb') as binfile:
        for i in range(genome_size):
            binfile.write((0).to_bytes(1, 'little'))
            if i % 1000000 == 0:
                print(f"{time.asctime()} | Processed {i/1000000} Mbases...", end="\r")

#------------------------------------------------------------------------------#
# Calculate intron, exon and promoter boundaries from the GFF file

if doGFF:
    print("# Parsing the GFF file...")
    with (gzip.open if genome_gff.endswith("gz") else open)(genome_gff, 'rt') as gff, \
                open(out_dir + "exons.coord", 'w') as exonfile, \
                open(out_dir + "introns.coord", 'w') as intronfile, \
                open(out_dir + "cisreg.coord", 'w') as cisregfile:

        line = next(gff)

        counter = 0

        while True:
            counter += 1
            if counter % 1000 == 0:
                print(f"{time.asctime()} | Processed {counter} annotations", end="\r")

            try:
                # If it's a gene, read every exon until the next gene:
                if "gene\t" in line:
                    row = line.strip("\n").split("\t")
                    gene_scaf, gene_start, gene_end, strand = row[0], int(row[3]), int(row[4]), row[6]
                    offset = abs_dict[gene_scaf]

                    if strand == "+":
                        promoter_start = offset + max(0, gene_start - 2000)
                        promoter_end = offset + gene_start

                    else:
                        promoter_start = offset + gene_end
                        promoter_end = offset + min(len_dict[gene_scaf], gene_end + 2000)

                    cisregfile.write("\t".join([str(promoter_start), str(promoter_end)]) + "\n")
                    #coordfile.write("\t".join(["gene", ex_scaf, str(gene_start), str(gene_end)]) + "\n")

                    exons = []
                    x = next(gff)
                    while "gene\t" not in x:
                        if "\texon\t" in x:
                            ex_row = x.strip("\n").split("\t")
                            exons.append([ex_row[0], int(ex_row[3]), int(ex_row[4])])
                        x = next(gff)

                    if len(exons) > 0:
                        if strand == "-":
                            exons.reverse()
                        for x in exons:
                            ex_scaf, ex_start, ex_end = x
                            if ex_start == gene_start and ex_end < gene_end:
                                exonfile.write("\t".join([str(offset + ex_start), str(offset + ex_end)]) + "\n")
                                intronfile.write("\t".join([str(offset + ex_end)]))
                            elif ex_start > gene_start and ex_end < gene_end:
                                intronfile.write("\t" + str(offset + ex_start) + "\n")
                                exonfile.write("\t".join([str(offset + ex_start), str(offset + ex_end)]) + "\n")
                                intronfile.write("\t".join([str(offset + ex_end)]))
                            elif ex_start > gene_start and ex_end == gene_end:
                                intronfile.write("\t" + str(offset + ex_start) + "\n")
                                exonfile.write("\t".join([str(offset + ex_start), str(offset + ex_end)]) + "\n")
                            elif ex_start == gene_start and ex_end == gene_end:
                                exonfile.write("\t".join([str(offset + ex_start), str(offset + ex_end)]) + "\n")
                        line = next(gff)

                else:
                    line = next(gff)

            except:
                break

    print("\n# Filling in exons...")
    feature_type = "exon"
    this_type = 2**codes[feature_type]
    with open(out_dir + "exons.coord", 'rt') as ifile, open(map_path, 'rb+') as binfile:
        for k, line in enumerate(ifile):
            start, end = [int(x) for x in line.strip("\n").split("\t")]
            length = end - start
            for i in range(length):
                binfile.seek(start + i)
                current = int.from_bytes(binfile.read(1), 'little')
                new = (current | this_type).to_bytes(1, 'little')
                binfile.seek(start + i)
                binfile.write(new)
            if k % 1000 == 0:
                print(f"{time.asctime()} | Processed {k} exons...", end="\r")

    print("\n# Filling in introns...")
    feature_type = "intron"
    this_type = 2**codes[feature_type]
    with open(out_dir + "introns.coord", 'rt') as ifile, open(map_path, 'rb+') as binfile:
        for k, line in enumerate(ifile):
            start, end = [int(x) for x in line.strip("\n").split("\t")]
            length = end - start
            for i in range(length):
                binfile.seek(start + i)
                current = int.from_bytes(binfile.read(1), 'little')
                new = (current | this_type).to_bytes(1, 'little')
                binfile.seek(start + i)
                binfile.write(new)
            if k % 1000 == 0:
                print(f"{time.asctime()} | Processed {k} introns...", end="\r")

    print("\n# Filling in cis-regulatory elements...")
    feature_type = "cisreg"
    this_type = 2**codes[feature_type]
    with open(out_dir + "cisreg.coord", 'rt') as ifile, open(map_path, 'rb+') as binfile:
        for k, line in enumerate(ifile):
            start, end = [int(x) for x in line.strip("\n").split("\t")]
            length = end - start
            for i in range(length):
                binfile.seek(start + i)
                current = int.from_bytes(binfile.read(1), 'little')
                new = (current | this_type).to_bytes(1, 'little')
                binfile.seek(start + i)
                binfile.write(new)
            if k % 1000 == 0:
                print(f"{time.asctime()} | Processed {k} cis-regulatory elements...", end="\r")

if doCPGI:
    print("\n# Filling in CpG clusters...")
    feature_type = "cpgi"
    this_type = 2**codes[feature_type]
    with (gzip.open if genome_cpgi.endswith("gz") else open)(genome_cpgi, 'rt') as ifile, open(map_path, 'rb+') as binfile:
        next(ifile)
        for k, line in enumerate(ifile):
            row = line.strip("\n").split("\t")
            scaf, start, end = row[0], int(row[1]), int(row[2])
            length = end - start
            offset = abs_dict[scaf]
            fullstart = offset + start
            for i in range(length):
                binfile.seek(fullstart + i)
                current = int.from_bytes(binfile.read(1), 'little')
                new = (current | this_type).to_bytes(1, 'little')
                binfile.seek(fullstart + i)
                binfile.write(new)
            if k % 1000 == 0:
                print(f"{time.asctime()} | Processed {k} CpG clusters...", end="\r")

if doCpG:
    print("\n# Filling in CpG sites...")
    feature_type = "cpg"
    this_type = 2**codes[feature_type]
    position = 0
    with (gzip.open if genome_file.endswith("gz") else open)(genome_file, 'rt') as ifile, open(map_path, 'rb+') as binfile:
        for line in ifile:
            if line.startswith(">"):
                pass
            else:
                for k, base in enumerate(line.strip("\n")):
                    if base.upper() == "C" and line[k+1].upper() == "G":
                        # Annotate the C
                        binfile.seek(position)
                        current = int.from_bytes(binfile.read(1), 'little')
                        new = (current | this_type).to_bytes(1, 'little')
                        binfile.seek(position)
                        binfile.write(new)
                        # Annotate the G
                        binfile.seek(position+1)
                        current = int.from_bytes(binfile.read(1), 'little')
                        new = (current | this_type).to_bytes(1, 'little')
                        binfile.seek(position+1)
                        binfile.write(new)
                    position += 1
            if position % 100000 == 0:
                print(f"{time.asctime()} | Processed {position/1000000} million base pairs...", end="\r")

#------------------------------------------------------------------------------#
# Annotating C/T SNPs from a VCF file

if doVCF:
    print("\n# Annotating C/T polymorphisms...")
    feature_type = "ct_snp"
    this_type = 2**codes[feature_type]

    i = 0
    with (gzip.open if dbsnp_path.endswith("gz") else open)(dbsnp_path, 'rt') as ifile, open(map_path, 'rb+') as binfile:
        for line in ifile:
            if not line.startswith("#"):
                row = line.split("\t")
                if all(row[x] in ("C", "T") for x in (3,4)) or all(row[x] in ("G", "A") for x in (3,4)):
                    scaf = row[0]
                    pos = int(row[1]) - 1
                    offset = abs_dict[scaf]
                    position = offset + pos
                    binfile.seek(position)
                    current = int.from_bytes(binfile.read(1), 'little')
                    new = (current | this_type).to_bytes(1, 'little')
                    binfile.seek(position)
                    binfile.write(new)
                    i += 1

                    if i > 0 and i % 100000 == 0:
                        print(f"{time.asctime()} | Processed {i/1000} thousand C/T polymorphisms...", end="\r")

#------------------------------------------------------------------------------#
# Writing a bin matrix with methylation data and feature flags

if doMatrix:

    print("\n# Preparing the binary methylation matrix...")

    with (gzip.open if matrix_file.endswith("gz") else open)(matrix_file, 'rt') as mfile, \
            open(map_path, 'rb') as bfile, \
            open(out_dir + "methylation_features.matrix", 'wb') as ofile:

        header = next(mfile).strip("\n").split("\t")
        base_header = [h.strip("_meth_cytosine") for h in header]
        sample_ids = [int((base_header.index(i)-1)/2) for i in samples]

        # Read in the first row:
        prev_row = next(mfile).strip("\n").split("\t")
        sp = prev_row[0].split(":")
        prev_fullpos = abs_dict[sp[0]] + int(sp[1]) - 1
        prev_count = [float(x) if x != "nan" else float(0) for x in prev_row[1:]]

        is_new_cpg = True
        # Start iterating over the file:
        for k, line in enumerate(mfile):
            this_row = line.strip("\n").split("\t")
            sp = this_row[0].split(":")
            this_fullpos = abs_dict[sp[0]] + int(sp[1]) - 1
            this_count = [float(x) if x != "nan" else float(0) for x in this_row[1:]]

            # if we are on the same scaffold, at 1bp distance:
            if this_fullpos == prev_fullpos + 1 and is_new_cpg:
                is_new_cpg = False
                data = []
                for i in sample_ids:
                    if not all(x == 0 for x in (prev_count[2 * i + 1], this_count[2 * i + 1])):
                        meth = (prev_count[2 * i] + this_count[2 * i]) / (prev_count[2 * i + 1] + this_count[2 * i + 1])
                    else:
                        meth = float(9.999999)
                    data.append(meth)

                bfile.seek(prev_fullpos)
                bflag_0 = bfile.read(1)
                bflag_1 = bfile.read(1)
                bflag = (int.from_bytes(bflag_0, 'little') | int.from_bytes(bflag_1, 'little')).to_bytes(1, 'little')
                s = struct.pack('f'*len(data), *data)
                ofile.write(bflag)
                ofile.write(s)
            else:
                is_new_cpg = True

            prev_row = next(mfile).strip("\n").split("\t")
            sp = prev_row[0].split(":")
            prev_fullpos = abs_dict[sp[0]] + int(sp[1]) - 1
            prev_count = [float(x) if x != "nan" else float(0) for x in prev_row[1:]]

        # for k, line in enumerate(mfile):
        #     row = line.strip("\n").split("\t")
        #     spos = row[0].split(":")
        #     scaf, pos = spos[0], int(spos[1])
        #     fullpos = abs_dict[scaf] + pos - 1
        #     data = [float(row[i]) for i in sample_ids]
        #     bfile.seek(fullpos)
        #     bflag = bfile.read(1)
        #     s = struct.pack('f'*len(data), *data)
        #     ofile.write(bflag)
        #     ofile.write(s)

            if k % 2000 == 0:
                print(f"{time.asctime()} | Processed {k/2} mCpG sites...", end="\r")

print("\n---------------------------------------------------------------------------------")
print(f"Completed in {str(datetime.timedelta(seconds=round(time.time() - start_time)))}")
