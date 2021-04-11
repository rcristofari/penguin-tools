import argparse

parser = argparse.ArgumentParser(description='Find sex-related scaffolds and identify sex of samples based on the distribution of heterozygosity. Input files are in the "012" format output by VCFtools. Windows are computed as numbers of SNPs, not as distances along the scaffolds, in order to keep weights homogeneous when computing the composite likelihoods.')
parser.add_argument('--geno', help='Basename of the input files (VCFtools option --012')
parser.add_argument('--out', help='Basename of the output files', default="Same as input prefix")
parser.add_argument('--wSize', help='Size of the SNP blocks used to compute heterozygosity [50]', default="50")
parser.add_argument('--wStep', help='Step between the SNP blocks [=wSize, non-overlapping]', default=None)
parser.add_argument('--alpha', help='p-value cutoff for bi-modality assignment [.001]', default=.001)
parser.add_argument('--sexes', help='Optional - A tab-separated files with sample names and sex assignments (homogametic / heterogametic)', default=None)
parser.add_argument('--chrom', help='Optional - Restrict processing to this chromosome / scaffold', default=None)
parser.add_argument('--doPlot', action='store_true', help='Output a plot of the classification, per scaffold')
parser.add_argument('--gff', help='Optional - a gff file for showing CDS intervals on the plots', default=None)
parser.add_argument('--verbose', action='store_true', help="Display all runtime / debugging details")
args = parser.parse_args()

# Required libraries
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import numpy.ma as ma
from sklearn.cluster import KMeans
from scipy.stats import norm
from scipy.stats.distributions import chi2

#########################################
# Load the individual and position files
posdf = pd.read_csv(args.geno + ".012.pos", delimiter = "\t", header=None, names=("chrom", "pos"))
ind = [x for x in pd.read_csv(args.geno + ".012.indv", delimiter = "\t", header=None)[0]]
nPos = posdf.shape[0]
nInd = len(ind)

# We extract the list of available scaffolds:
scaffolds = [x for x in sorted(set(posdf["chrom"]))]

#########################################
# Set the window size and step:
wSize = int(args.wSize)
if not args.wStep:
    wStep = int(args.wSize)

#########################################
# Load the genotype file:
genoList = []
with open(args.geno + ".012") as ifile:
    for line in ifile:
        genoList.append([int(x) for x in line.strip("\n").split("\t")])
geno = np.array(genoList)

#########################################
# Load the gff file:
if args.gff:
    gff = pd.read_csv(args.gff, names=('scaffold', 'source', "feature", "start", "end", "score", "strand", "phase", "attributes"), delimiter="\t")
    gff = gff.loc[gff["feature"] == "mRNA"]
else:
    gff = None

#########################################
# Load the sex assignment file:
if args.sexes:
    sexes = pd.read_csv(args.sexes, delimiter="\t", header=None, names=("ind", "sex"))
    sexes = dict(zip(sexes["ind"], sexes["sex"]))
    sexAssign = [0 if sexes[x] == "homogametic" else 1 for x in ind]
else:
    sexes, sexAssign = None, None

#########################################
# A function to distinguish uni- vs bi-modal distributions using K-means clustering and a likelihood ratio test
def k2_lrt(dat, init_="k-means++", alpha=float(args.alpha)):
    # init_ is a vector of 0 and 1, describing the a priori groups for initialisation

    # Perform K=2 clustering:
    dat = np.array(dat)
    dat = dat.reshape(-1, 1)

    # Parse the initialisation parameters
    if type(init_) is list and len(init_) == len(dat):
        init_ = np.array(init_)
        mean_0 = np.mean(dat[np.where(init_ == 0)])
        mean_1 = np.mean(dat[np.where(init_ == 1)])
        init_ = np.array([[mean_0], [mean_1]])
        kmeans = KMeans(init=init_, n_clusters=2, n_init=1, max_iter=300, random_state=42)
    else:
        kmeans = KMeans(init="k-means++", n_clusters=2, n_init=10, max_iter=300, random_state=42)

    kmeans.fit(dat)
    # One cluster model log-likelihood:
    mu_1 = np.mean(dat)
    sigma_1 = np.std(dat)
    lnL_1 = np.sum(np.log(norm.pdf(dat, mu_1, sigma_1)))

    # Two cluster model log-likelihood:
    # Cluster 0:
    dat_0 = dat[np.where(kmeans.labels_ == 0)]
    mu_2_0 = np.mean(dat_0)
    sigma_2_0 = np.std(dat_0)
    lnL_2_0 = np.sum(np.log(norm.pdf(dat_0, mu_2_0, sigma_2_0)))

    # Cluster 1:
    dat_1 = dat[np.where(kmeans.labels_ == 1)]
    mu_2_1 = np.mean(dat_1)
    sigma_2_1 = np.std(dat_1)
    lnL_2_1 = np.sum(np.log(norm.pdf(dat_1, mu_2_1, sigma_2_1)))

    # Likelihood-ratio test:
    lnL_2 = lnL_2_0 + lnL_2_1
    LRT = - 2 * (lnL_1 - lnL_2)
    pval = chi2.sf(LRT, 2 + len(dat))
    # Degrees of freedom: I consider the base model has 1 df (the mean)

    # If we declare the test significant:
    if pval <= alpha:
        isBimodal = 1
        # We polarise the clusters: 0 for low het, 1 for high het
        if mu_2_0 >= mu_2_1:
            assign = [1 if x == 0 else 0 for x in kmeans.labels_]
        else:
            assign = [x for x in kmeans.labels_]
    else:
        isBimodal = 0
        assign = [0] * dat.shape[0]

    return ([assign, pval, isBimodal])

#########################################
# Compute heterozygosity on blocks of SNPs:
def windowHet(geno, posdf, scaf, wSize=50, wStep=None):

    if not wStep:
        wStep = wSize

    these_idx = [x for x in posdf.loc[posdf["chrom"] == scaf].index]
    these_pos = np.array([x for x in posdf["pos"].loc[these_idx]])
    these_geno = geno[:, these_idx]

    # Bases with "0" and "2" are equally monomorphic. Here we set them all to "0":
    these_geno[these_geno == 2] = 0

    nGeno = these_geno.shape[1]
    nWindows = ((nGeno - wSize) // wStep) + 2

    these_means = np.full([nInd, nWindows], np.nan)
    wCenter = []
    wCenterPos = []
    start = 0
    stop = wSize
    w = 0

    while stop <= (nGeno + wStep):
        wData = these_geno[:, start:stop]
        masked = ma.array(wData, mask=[wData == -1])
        these_means[:, w] = masked.mean(axis=1)
        wCenter.append(np.mean([start, stop]))
        if stop >= nGeno:
            realStop = nGeno - 1
        else:
            realStop = stop
        wCenterPos.append(np.mean([these_pos[start], these_pos[realStop]]))

        start += wStep
        stop += wStep
        w += 1

    return ([these_means, wCenterPos])

#########################################
# Classify blocks in uni- or bimodal:
def classifyWindows(meanArray, centers, init=None, gff=None, scaf=None, doPlot=True):
    assignArray = np.full([meanArray.shape[0], meanArray.shape[1]], np.nan)
    pvals, whichBimodals = [], []

    if doPlot:
        plt.rcParams['figure.figsize'] = [18, 6]

        if gff is not None and scaf is not None:
            # Draw mRNA bands:
            this_gff = gff.loc[gff["scaffold"] == scaf]

            for i in range(this_gff.shape[0]):
                plt.axvspan(this_gff.iloc[i]["start"], this_gff.iloc[i]["end"], color="k", alpha=.1)

    for x in range(meanArray.shape[1]):
        assign, pval, isBimodal = k2_lrt(meanArray[:, x], init_=init, alpha=.01)
        assignArray[:, x] = assign
        pvals.append(pval)
        whichBimodals.append(isBimodal)

        if doPlot:
            if isBimodal == 1:
                col = ["#61E294" if a == 1 else "#FFC857" for a in assign]
            else:
                col = "#62929E"
            plt.scatter([centers[x]] * meanArray.shape[0], meanArray[:, x], c=col, s=2)

    if doPlot:
        plt.show()

    return ([assignArray, pvals, whichBimodals])

#########################################
# Calculate the likelihood for being a Z-linked scaffold:
def ScafLRT(groups, meanArray):
    # Note that this will ONLY be accurate with NON-OVERLAPPING windows..!
    from scipy.stats import norm
    from scipy.stats.distributions import chi2

    K1_lnL = []
    K2_lnL = []

    for m in range(meanArray.shape[1]):
        dat = meanArray[:, m]

        # One cluster model log-likelihood:
        mu_1 = np.mean(dat)
        sigma_1 = np.std(dat)
        lnL_1 = np.sum(np.log(norm.pdf(dat, mu_1, sigma_1)))
        K1_lnL.append(lnL_1)

        # Two cluster model log-likelihood:
        # Cluster 0:
        dat_0 = dat[np.where(np.array(groups) == 0)]
        mu_2_0 = np.mean(dat_0)
        sigma_2_0 = np.std(dat_0)
        lnL_2_0 = np.sum(np.log(norm.pdf(dat_0, mu_2_0, sigma_2_0)))

        # Cluster 1:
        dat_1 = dat[np.where(np.array(groups) == 1)]
        mu_2_1 = np.mean(dat_1)
        sigma_2_1 = np.std(dat_1)
        lnL_2_1 = np.sum(np.log(norm.pdf(dat_1, mu_2_1, sigma_2_1)))
        lnL_2 = lnL_2_0 + lnL_2_1

        K2_lnL.append(lnL_2)

    # Likelihood-ratio test:
    LRT = - 2 * (np.sum(K1_lnL) - np.sum(K2_lnL))
    pval = chi2.sf(LRT, 2 * meanArray.shape[1] - 2 + len(dat))
    # One degree of freedom for each additional mean and sd, and one per a priori individual assignment

    return (pval)

# If we process only one scaffold:
if args.chrom:
    meanArray, centers = windowHet(geno, posdf, scaf=args.chrom, wSize=wSize, wStep=wStep)

    if not args.sexes:
        init = "k-means++"
    else:
        init = sexAssign
    assignArray, pvals, whichBimodals = classifyWindows(meanArray=meanArray, centers=centers, init=init, gff=gff, doPlot=args.doPlot)

    if not args.sexes:
        groups = [int(round(x)) for x in np.mean(assignArray, axis=1)]
    else:
        groups = sexAssign
    scaf_pval = ScafLRT(groups, meanArray)

    print(scaf_pval)


# TODO:
# Test how well kmeans converge with a false init assignment
# Loop through scaffolds, and use the first significantly sex-linked scaffold for initial sex assignment
# Reassess the assignment against the next significant sex-linked scaffold. This will work if kmeans is robuts to init state.
# Return pvalue for the proposed grouping, but also pvalue for an eventual better grouping:
# if assignments are provided, also compute the scaffold likelihood based on the mean assignments and remember both
# Return a table of scaffolds and pvalues, as well as a list of sex assignments, and Z-linked chromosomes

# At the scaffold level, Benjamini-Hochberg correction for non-independent multiple tests