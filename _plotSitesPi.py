# This reads in scaffold stats from a vcf file (gzipped or not) and the sites pi file from vcftools, and produces a plot
# Use: plotSitesPi("/path/to/vcf.gz", "path/to/file.sites.pi")

def plotSitesPi(vcf_path, sitesPi_path):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import gzip

    # Load the scaffold size dictionary from the VCF file:
    scaf, length = [], []
    with (gzip.open if vcf_path.endswith(".gz") else open)(vcf_path, "rt") as ifile:
        for line in ifile:
            if "##contig=" in line:
                row = line.strip("\n").split("=")
                length.append(int(row[3].strip(">")))
                scaf.append(row[2].split(",")[0])
            elif "#CHROM" in line:
                break
    scafDict = dict(zip(scaf, length))


    pi = pd.read_csv(sitesPi_path, delimiter="\t")

    # Sort the dataframe by scaffold length:
    pi["LEN"] = [scafDict[s] for s in pi["CHROM"]]
    pi = pi.sort_values(['LEN', 'CHROM', 'POS'], ascending=(False, True, True))
    pi.reset_index(inplace=True, drop=True)

    # Add a FULLPOS column corresponding to position + cumulative sum of sorted scaffold lengths
    uniquePairs = [x for x in sorted(set([x for x in zip(pi["CHROM"], pi["LEN"])]), key=lambda t: t[1], reverse=True)]
    s = [x[0] for x in uniquePairs]
    l = [x[1] for x in uniquePairs]
    l = [0] + l[:-1]
    csDict = dict(zip(s, np.cumsum(l)))

    pi["FULLPOS"] = [csDict[x] + pi["POS"][i] for i, x in enumerate(pi["CHROM"])]

    # Draw the boundaries of bands for scaffolds:
    stop = []
    for i, c in enumerate(pi["CHROM"][:-1]):
        if c != pi["CHROM"][i + 1]:
            stop.append(i)

    start = [x + 1 for x in stop[:-1]]
    start = [0] + start

    # for color bands, we only want every second band:
    bands = [x for i, x in enumerate(zip(start, stop)) if (i % 2 == 0)]

    # Redraw bands by position and not index:
    bandsPOS = []
    for b in bands:
        b0 = pi["FULLPOS"][b[0]]
        b1 = pi["FULLPOS"][b[1]]
        bandsPOS.append((b0, b1))

    # Get scaffold median:
    pi_median, pi_loc = [], []
    for i, s in enumerate(start):
        pi_median.append(np.nanmedian(pi["PI"][s:stop[i]]))
        pi_loc.append(np.mean([pi["FULLPOS"][s], pi["FULLPOS"][stop[i]]]))

    plt.rcParams['figure.figsize'] = [20, 6]
    for b in bandsPOS:
        plt.axvspan(b[0], b[1], facecolor='g', alpha=0.25)
    plt.scatter(pi["FULLPOS"], pi["PI"], s=.1)
    plt.scatter(pi_loc, pi_median, c="r", marker="+", s=10)
    plt.plot(pi_loc, pi_median, c="r", lw=.5)
    plt.show()
