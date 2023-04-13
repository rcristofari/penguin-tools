import argparse

parser = argparse.ArgumentParser(description="Create a feature map for a reference genome from a GFF file and a CpG cluster map. Create a binary methylation matrix file referenced to genome features, for a given set of individuals.")
parser.add_argument('--ref', help='Reference genome path in (gz)fasta')
parser.add_argument('--gff', help='Annotation file path (gz)gff')
parser.add_argument('--cpgi', help='CpG cluster file (gz)txt file from gCluster')
parser.add_argument('--vcf', help='SNP database in (gz)VCF format, to annotate C/T polymorphisms')
parser.add_argument('--cpg', help='Annotate CpG sites from the reference genome', action='store_true')
parser.add_argument('--map', help='Binary feature map file output by a previous run')
parser.add_argument('--matrix', help='Methylation count matrix from bsbolt, mC and totalC per sample')
parser.add_argument('--samples', help='Sample list (text file, one sample per row)')
parser.add_argument('--output', help='Output folder path')
parser.add_argument('--calls', help='Output individual methylation call files', action='store_true')
parser.add_argument('--suffix', help="Suffix to identify output files")
parser.add_argument('--includes', help="Features to be included in the output files (as integer)")
parser.add_argument('--excludes', help="Features to be excluded from the output files (as integer)")
parser.add_argument('--mincov', help="Minimum avg coverage at the CpG site to include a sample-site [0]")
parser.add_argument('--maxcov', help="Maximum avg coverage at the CpG site to include a sample-site [inf]")
parser.add_argument('--missing', help="Maximum proportion of missing individuals after coverage filtering [1]")
parser.add_argument('--stdev', help="Minimum standard deviation at the filtered site to keep it in output [0]")
parser.add_argument('--skip-single', help="Skip sites where only one base is covered (e.g. non-CpG sites, or RRBS data)", action='store_true')
parser.add_argument('--do-glm', help="Output a file ready for loading as an R GLM input", action = 'store_true')
parser.add_argument('--traits', help="A file containing individual traits, to add in the GLM input file")
args = parser.parse_args()

import struct, os, gzip, time, datetime
import numpy as np

start_time = time.time()

print("---------------------------------------------------------------------------------")
#----------------------------------------------------------------------------------------------------------------------#
# The exponents used for the bit codes:
codes = {'cpg':0,
         'cpgi':1, #2
         'cisreg':2, #4
         'intron':3, #8
         'exon_first':4, #16
         'exon_other':5, #32
         'ct_snp':6} #64
# ADD CpG SHORES AND SHELVES !

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
if args.gff:
    doGFF = True
    genome_gff = args.gff
    if not os.path.isfile(genome_gff):
        print(f"GFF file {genome_gff} could not be found")
        quit()
else:
    doGFF = False

#----------------------------------------------------------------------------------------------------------------------#
if args.cpgi:
    doCPGI = True
    genome_cpgi = args.cpgi
    if not os.path.isfile(genome_cpgi):
        print(f"CpG cluster file {genome_cpgi} could not be found")
        quit()
else:
    doCPGI = False

#----------------------------------------------------------------------------------------------------------------------#
if args.vcf:
    doVCF = True
    dbsnp_path = args.vcf
    if not os.path.isfile(dbsnp_path):
        print(f"VCF SNP file {dbsnp_path} could not be found")
        quit()
else:
    doVCF = False

#----------------------------------------------------------------------------------------------------------------------#
if args.cpg:
    doCpG = True
else:
    doCpG = False

#----------------------------------------------------------------------------------------------------------------------#
if args.output:
    out_dir = args.output
    if not os.path.isdir(out_dir):
        print(f"Output directory {out_dir} could not be found, defaulting to current directory")
        out_dir = "./"
    if not out_dir.endswith("/"):
        out_dir += "/"
else:
    out_dir = "./"
    print(f"Defaulting to current directory {os.path.abspath(os.getcwd())}")

#----------------------------------------------------------------------------------------------------------------------#
if args.map:
    doMap = False
    map_path = args.map
    if not os.path.isfile(map_path):
        print(f"Map file {map_path} could not be found")
        quit()
else:
    doMap = True
    map_path = out_dir + "features.map"

#----------------------------------------------------------------------------------------------------------------------#
if args.matrix:
    doMatrix = True
    matrix_file = args.matrix
    if not os.path.isfile(matrix_file):
        print(f"Matrix file {matrix_file} could not be found")
        quit()

    if args.samples:
        sample_file = args.samples
        if not os.path.isfile(sample_file):
            print(f"Sample list {sample_file} could not be found")
            quit()
        else:
            samples = []
            with open(sample_file, 'r') as ifile:
                for line in ifile:
                    if not line.startswith("#"):
                        samples.append(line.strip("\n"))
            print(f"Processing data for {len(samples)} samples.")
    else:
        # Get all samples from the header of the methylation matrix:
        with (gzip.open if matrix_file.endswith("gz") else open)(matrix_file, 'rt') as mfile:
            line = next(mfile).strip("\n").split("\t")
            samples = [x.strip("_meth_cytosine") for x in line if x.endswith("_meth_cytosine")]
            print(f"Processing data for all {len(samples)} samples.")
else:
    doMatrix = False

#----------------------------------------------------------------------------------------------------------------------#
if args.mincov:
    try:
        mincov = int(args.mincov)
        print(f"Filtering out genotypes with a coverage <= {mincov} X")
    except ValueError:
        mincov = None
        print(f"Invalid input {args.mincov} for minimum coverage")
        quit()
else:
    mincov = 0

#----------------------------------------------------------------------------------------------------------------------#
if args.maxcov:
    try:
        maxcov = int(args.maxcov)
        print(f"Filtering out genotypes with a coverage >= {maxcov} X")
    except ValueError:
        maxcov = None
        print(f"Invalid input {args.maxcov} for maximum coverage")
        quit()
else:
    maxcov = float('inf')

#----------------------------------------------------------------------------------------------------------------------#
if args.missing:
    try:
        missing = int(round(float(args.missing)*len(samples)))
        print(f"Removing sites missing in at least {int(round(float(args.missing)*100))}% of samples ({missing} samples)")
    except ValueError:
        missing = None
        print(f"Invalid input {args.missing} for maximum coverage")
        quit()
    if float(args.missing) < 0 or float(args.missing) > 1:
        print(f"Proportion of missing samples must be between 0 and 1 ({args.missing})")
        quit()
else:
    missing = len(samples)

#----------------------------------------------------------------------------------------------------------------------#
if args.stdev:
    try:
        stdev = float(args.stdev)
        print(f"Removing sites with a standard deviation below {stdev}")
    except ValueError:
        stdev = None
        print(f"Invalid input {args.stdev} for minimum standard deviation")
        quit()
else:
    stdev = None

#----------------------------------------------------------------------------------------------------------------------#
if args.calls:
    doCalls = True
    print("Making individual call files")
else:
    doCalls = False

#----------------------------------------------------------------------------------------------------------------------#
if args.do_glm:
    doGLM = True
    print("Outputting GLM input file")
else:
    doGLM = False

#----------------------------------------------------------------------------------------------------------------------#
if args.traits:
    traits_file = args.traits
    if not os.path.isfile(traits_file):
        print(f"Traits file {traits_file} could not be found")
        quit()
    else:
        print(f"Individual traits file: {args.traits}")
else:
    traits_file = None

#----------------------------------------------------------------------------------------------------------------------#
if args.includes:
    includes = int(args.includes)
else:
    includes = 0

#----------------------------------------------------------------------------------------------------------------------#
if args.excludes:
    excludes = int(args.excludes)
else:
    excludes = 0

#----------------------------------------------------------------------------------------------------------------------#
if args.skip_single:
    skip_single = True
    print("Skipping Cs from incomplete CpG or non-CpG sites")
else:
    skip_single = False

#----------------------------------------------------------------------------------------------------------------------#
if args.suffix:
    suffix = args.suffix
else:
    if includes + excludes > 0:
        suffix = f"{includes}_not_{excludes}"
    else:
        suffix = "allSites"

#----------------------------------------------------------------------------------------------------------------------#
if includes > 0:
    include_str = [c for c in codes if (1<<codes[c] & includes) > 0]
    print("Including " + ", ".join(include_str))
if excludes > 0:
    exclude_str = [c for c in codes if (1<<codes[c] & excludes) > 0]
    print("Excluding " + ", ".join(exclude_str))

print("---------------------------------------------------------------------------------")

#----------------------------------------------------------------------------------------------------------------------#
# Figure out the total length of the reference genome:

print("# Calculating target genome size...", end='')
genome_size = 0
with (gzip.open if genome_file.endswith("gz") else open)(genome_file, 'rt') as gfile:
    for line in gfile:
        if line.startswith(">"):
            pass
        else:
            genome_size += len(line.strip('\n'))
print(f"done. Genome size = {genome_size} bp")

#------------------------------------------------------------------------------#
# Create a coordinate reference file for the scaffolds in the genome:

# This is zero-based

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

#----------------------------------------------------------------------------------------------------------------------#
# Create the feature map file

if doMap:
    print("# Initializing the feature map file...", end="\r")
    # Initialise an empty code file:
    with open(map_path, 'wb') as binfile:
        binfile.write((0).to_bytes(genome_size, 'little'))
        print("done.")

#----------------------------------------------------------------------------------------------------------------------#
# Calculate intron, exon and promoter boundaries from the GFF file
# GFF file is zero-based

if doGFF:
    print("# Parsing the GFF file...")
    with (gzip.open if genome_gff.endswith("gz") else open)(genome_gff, 'rt') as gff, \
                open(out_dir + "exons_first.coord", 'w') as first_exonfile, \
                open(out_dir + "exons_other.coord", 'w') as other_exonfile, \
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
                        for e, x in enumerate(exons):
                            ex_scaf, ex_start, ex_end = x

                            # Point to the right output file to separate the first exon from the others:
                            if e == 0:
                                exonfile = first_exonfile
                            else:
                                exonfile = other_exonfile

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

    print("\n# Filling in first exons...")
    feature_type = "exon_first"
    this_type = 2**codes[feature_type]
    with open(out_dir + "exons_first.coord", 'rt') as ifile, open(map_path, 'rb+') as binfile:
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
                print(f"{time.asctime()} | Processed {k} exons", end="\r")

    print("\n# Filling in other exons...")
    feature_type = "exon_other"
    this_type = 2**codes[feature_type]
    with open(out_dir + "exons_other.coord", 'rt') as ifile, open(map_path, 'rb+') as binfile:
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
                print(f"{time.asctime()} | Processed {k} exons", end="\r")

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
                print(f"{time.asctime()} | Processed {k} introns", end="\r")

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
                print(f"{time.asctime()} | Processed {k} cis-regulatory elements", end="\r")

#----------------------------------------------------------------------------------------------------------------------#
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
                print(f"{time.asctime()} | Processed {k} CpG clusters", end="\r")

#----------------------------------------------------------------------------------------------------------------------#
if doCpG:
    print("\n# Filling in CpG sites...")
    feature_type = "cpg"
    this_type = 2**codes[feature_type]
    position = 0
    with (gzip.open if genome_file.endswith("gz") else open)(genome_file, 'rt') as ifile, open(map_path, 'rb+') as binfile:

        for line in ifile:
            if line.startswith(">"):
                last_is_C = False

            else:
                # Handle CpG overhanging 2 lines
                if line[0].upper() == "G" and last_is_C:
                    # Annotate the C
                    binfile.seek(position-1) #
                    current = int.from_bytes(binfile.read(1), 'little')
                    new = (current | this_type).to_bytes(1, 'little')
                    binfile.seek(position-1) #
                    binfile.write(new)
                    # Annotate the G
                    binfile.seek(position) #
                    current = int.from_bytes(binfile.read(1), 'little')
                    new = (current | this_type).to_bytes(1, 'little')
                    binfile.seek(position) #
                    binfile.write(new)

                if line[-2].upper() == "C":
                    last_is_C = True
                else:
                    last_is_C = False

                for k, base in enumerate(line[:-1]):
                    if base.upper() == "C" and line[k + 1].upper() == "G":
                        # Annotate the C
                        binfile.seek(position)
                        current = int.from_bytes(binfile.read(1), 'little')
                        new = (current | this_type).to_bytes(1, 'little')
                        binfile.seek(position)
                        binfile.write(new)
                        # Annotate the G
                        binfile.seek(position + 1)
                        current = int.from_bytes(binfile.read(1), 'little')
                        new = (current | this_type).to_bytes(1, 'little')
                        binfile.seek(position + 1)
                        binfile.write(new)
                    position += 1

                    if position % 100000 == 0:
                        print(f"{time.asctime()} | Processed {position / 1000000} million base pairs", end="\r")

#----------------------------------------------------------------------------------------------------------------------#
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
                        print(f"{time.asctime()} | Processed {i/1000} thousand C/T polymorphisms", end="\r")

#----------------------------------------------------------------------------------------------------------------------#
# Writing a bin matrix with methylation data and feature flags

if doMatrix:

    print("\n# Preparing the binary methylation matrix...")

    with (gzip.open if matrix_file.endswith("gz") else open)(matrix_file, 'rt') as mfile, \
            open(map_path, 'rb') as bfile, \
            open(out_dir + f"methylation_features.{suffix}.matrix", 'wb') as ofile:

        # Counter to display the number of retained sites
        kept, kept_single, lowstd = 0, 0, 0

        header = next(mfile).strip("\n").split("\t")
        base_header = [h.strip("_meth_cytosine") for h in header]
        sample_ids = [int((base_header.index(i)-1)/2) for i in samples]

        # If we are outputting individual call files, initialize these:
        if doCalls:
            cfiles = []
            for s in samples:
                cfiles.append(open(out_dir + s + f".{suffix}.calls", 'w'))
            for f in cfiles:
                f.write("chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n")

        # If we are outputting a GLM file, initialize it:
        if doGLM:
            glmfile = open(out_dir + f"methylation.{suffix}.glm", 'w')

            # Prepare a dict of dicts for the traits file, to enrich the GLM file
            if traits_file:
                with open(traits_file, 'r') as tfile:
                    header = next(tfile).strip("\n").split("\t")
                    ntraits = len(header) - 1
                    trait_names = header[1:]
                    trait_list = [[] for x in range(ntraits)]
                    trait_samples = []
                    for line in tfile:
                        row = line.strip("\n").split("\t")
                        trait_samples.append(row[0])
                        for x, l in enumerate(row[1:]):
                            trait_list[x].append(l)
                    trait_dicts = [[] for x in range(ntraits)]

                    for x, l in enumerate(trait_list):
                        trait_dicts[x] = dict(zip(trait_samples, l))
                    traits = dict(zip(trait_names, trait_dicts))

            if traits_file:
                glmfile.write("chrom\tpos\tfullpos\tsample\tC\tT\tfeature\tcpgi")
                for t in trait_names:
                    glmfile.write("\t" + t)
                glmfile.write("\n")
            else:
                glmfile.write("chrom\tpos\tfullpos\tsample\tC\tT\tfeature\tcpgi\n")

        # Read in the first row:
        prev_row = next(mfile).strip("\n").split("\t")
        prev_scaf, prev_pos = prev_row[0].split(":")
        prev_fullpos = abs_dict[prev_scaf] + int(prev_pos) - 1
        prev_count = [float(x) if x != "nan" else float(0) for x in prev_row[1:]]
        is_new_cpg = False
        doWrite = False

        # Start iterating over the file:
        for k, line in enumerate(mfile):
            this_row = line.strip("\n").split("\t")
            this_scaf, this_pos = this_row[0].split(":")
            this_fullpos = abs_dict[this_scaf] + int(this_pos) - 1
            this_count = [float(x) if x != "nan" else float(0) for x in this_row[1:]]

            # CASE OF THE SECOND BASE IN A CpG SITE: WE AVERAGE AND WRITE BOTH TO 1ST POSITION
            if this_scaf == prev_scaf and this_fullpos == prev_fullpos + 1 and not is_new_cpg:

                is_new_cpg = True

                data, coverage, CT = [], [], []
                for i in sample_ids:
                    if not all(x == 0 for x in (prev_count[2 * i + 1], this_count[2 * i + 1])):
                        meth = (prev_count[2 * i] + this_count[2 * i]) / (prev_count[2 * i + 1] + this_count[2 * i + 1])
                        cov = prev_count[2 * i + 1] + this_count[2 * i + 1]
                        ct = [(prev_count[2 * i] + this_count[2 * i]), (cov - (prev_count[2 * i] + this_count[2 * i]))]
                    else:
                        meth = float(9.999999)
                        cov = 0
                        ct = [0, 0]

                    data.append(meth)
                    coverage.append(cov)
                    CT.append(ct)

                bfile.seek(prev_fullpos)
                bflag_0 = bfile.read(1)
                bflag_1 = bfile.read(1)
                bpos = (prev_fullpos).to_bytes(4, 'little')
                flag = int.from_bytes(bflag_0, 'little') | int.from_bytes(bflag_1, 'little')

                if ((flag & includes) == includes) and ((flag & excludes) == 0):

                    # Here, do the filtering:
                    masked = [d if (coverage[i]>=2*mincov and coverage[i]<=2*maxcov) else float(9.999999) for i, d in enumerate(data)]

                    if sum(m == float(9.999999) for m in masked) <= missing:
                        kept += 1
                        bflag = (flag).to_bytes(1, 'little')
                        doWrite = True

                        if stdev:
                            m_array = np.array(masked)
                            m_array[np.where(m_array > 1)] = np.nan
                            this_stdev = np.nanstd(m_array)
                            if this_stdev < stdev:
                                doWrite = False
                                kept -= 1
                                lowstd += 1

            # CASE WHERE WE HAVE AN ISOLATED mC SITE: WE WRITE IT OUT IF NOT DISABLED
            elif this_scaf == prev_scaf and not is_new_cpg and not skip_single:

                data, coverage, CT = [], [], []
                for i in sample_ids:
                    if prev_count[2 * i + 1] != 0:
                        meth = (prev_count[2 * i]) / (prev_count[2 * i + 1])
                        cov = prev_count[2 * i + 1]
                        ct = [prev_count[2 * i], (cov - prev_count[2 * i])]
                    else:
                        meth = float(9.999999)
                        cov = 0
                        ct = [0, 0]

                    data.append(meth)
                    coverage.append(cov)
                    CT.append(ct)

                bfile.seek(prev_fullpos)
                bflag = bfile.read(1)
                bpos = (prev_fullpos).to_bytes(4, 'little')
                flag = int.from_bytes(bflag, 'little')

                if ((flag & includes) == includes) and ((flag & excludes) == 0):

                    # Here, do the filtering:
                    masked = [d if (coverage[i]>=2*mincov and coverage[i]<=2*maxcov) else float(9.999999) for i, d in enumerate(data)]

                    if sum(m == float(9.999999) for m in masked) <= missing:
                        kept_single += 1
                        doWrite = True

                        if stdev:
                            m_array = np.array(masked)
                            m_array[np.where(m_array > 1)] = np.nan
                            this_stdev = np.nanstd(m_array)
                            if this_stdev < stdev:
                                doWrite = False
                                kept_single -= 1
                                lowstd += 1

            else:
                is_new_cpg = False

            if doWrite:

                s = struct.pack('f' * len(masked), *masked)
                ofile.write(bpos)
                ofile.write(bflag)
                ofile.write(s)

                # If we are writing call files, write the data to each file:
                if doCalls:
                    for c, f in enumerate(cfiles):
                        schrom, spos = prev_row[0].split(":")
                        percent0 = (round(masked[c] * 100, 2) if masked[c] != float(9.999999) else 0)
                        percent1 = (round((1-masked[c]) * 100, 2) if masked[c] != float(9.999999) else 0)
                        f.write("\t".join([prev_row[0], schrom, spos, "F", str(int(coverage[c])), str(percent0), str(percent1)]) + "\n")

                if doGLM:
                    site_features = [c for c in codes if (1 << codes[c] & flag) > 0]
                    if 'cpgi' in site_features:
                        is_cpgi = 'Y'
                    else:
                        is_cpgi = 'N'
                    if 'cisreg' in site_features:
                        this_feature = 'cisreg'
                    elif 'exon_first' in site_features:
                        this_feature = 'exon_first'
                    elif 'exon_other' in site_features:
                        this_feature = 'exon_other'
                    elif 'intron' in site_features:
                        this_feature = 'intron'
                    else:
                        this_feature = 'intergenic'

                    for y, z in enumerate(CT):
                        glmline = [str(prev_scaf),
                                   str(prev_pos),
                                   str(prev_fullpos),
                                   str(samples[y])] + [str(int(b)) for b in CT[y]] + [this_feature,
                                   is_cpgi]
                        if traits_file:
                            for t in traits:
                                glmline += [str(traits[t][samples[y]])]
                        glmfile.write("\t".join(glmline) + "\n")

                doWrite = False



            prev_row, prev_fullpos, prev_count, prev_scaf, prev_pos = this_row, this_fullpos, this_count, this_scaf, this_pos

            if k % 2000 == 0:
                retaining_rate = round(((int(kept*2 + kept_single))/int(k))*100, 2) if int(k) > 0 else 0
                print(f"{time.asctime()} | Processed {int(k)} sites | kept {int(kept*2 + kept_single)} bases | retained {int(kept)} CpG sites and {int(kept_single)} single sites | removed {lowstd} low-variance sites | Retaining rate {retaining_rate}%", end="\r")

    if doCalls:
        for f in cfiles:
            f.close()

    if doGLM:
        glmfile.close()

print("\n---------------------------------------------------------------------------------")
print(f"Completed in {str(datetime.timedelta(seconds=round(time.time() - start_time)))}")
