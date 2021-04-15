import argparse

parser = argparse.ArgumentParser(description='Find sex-related scaffolds and identify sex of samples based on the distribution of heterozygosity. Input files are in the "012" format output by VCFtools. Windows are computed as numbers of SNPs, not as distances along the scaffolds, in order to keep weights homogeneous when computing the composite likelihoods.')
parser.add_argument('--geno', help='Basename of the input files (VCFtools option --012')
parser.add_argument('--out', help='Basename of the output files', default="Same as input prefix")
parser.add_argument('--wSize', help='Size of the SNP blocks used to compute heterozygosity [50]', default="50")
parser.add_argument('--wStep', help='Step between the SNP blocks [=wSize, non-overlapping]', default=None)
parser.add_argument('--minSNPs', help='Minimum number of SNPs to process the scaffold [10]', default=10)
parser.add_argument('--alpha', help='p-value cutoff for bi-modality assignment [.001]', default=.001)
parser.add_argument('--sexes', help='Optional - A tab-separated files with sample names and sex assignments (homogametic / heterogametic)', default=None)
parser.add_argument('--chrom', help='Optional - Restrict processing to this chromosome / scaffold', default=None)
parser.add_argument('--doPlot', action='store_true', help='Output a plot of the classification, per scaffold')
parser.add_argument('--gff', help='Optional - a gff file for showing CDS intervals on the plots', default=None)
parser.add_argument('--verbose', action='store_true', help="Display all runtime / debugging details")
args = parser.parse_args()

# Required libraries
import os, math
import pandas as pd
import numpy as np
import numpy.ma as ma
from collections import Counter
from sklearn.cluster import KMeans
from scipy.stats import norm
from scipy.stats.distributions import chi2

if args.doPlot:
    import matplotlib.pyplot as plt

###########################################################
## REMOVE THIS LATER AND CATCH THE SPECIFIC WARNING !!!!  #
import warnings                                           #
warnings.filterwarnings("ignore")                         #
###########################################################

#########################################
# Load the individual and position files
print("------------------------------------------------")
print("Loading SNP metadata from " + os.path.abspath(args.geno + ".012.pos"))
posdf = pd.read_csv(args.geno + ".012.pos", delimiter = "\t", header=None, names=("chrom", "pos"))
ind = [x for x in pd.read_csv(args.geno + ".012.indv", delimiter = "\t", header=None)[0]]
nPos = posdf.shape[0]
nInd = len(ind)
print("Read " + str(posdf.shape[0]) + " positions")
print("------------------------------------------------")

#########################################
# We extract the list of available scaffolds:
nSNPs = Counter(posdf["chrom"])
scaffolds = sorted([s for s in nSNPs if nSNPs[s] > int(args.minSNPs)])


#########################################
# Set the window size and step:
wSize = int(args.wSize)
if not args.wStep:
    wStep = int(args.wSize)
else:
    wStep = args.wStep

#########################################
# Load the genotype file:
print("Loading genotypes from " + os.path.abspath(args.geno + ".012"))
genoList = []
with open(args.geno + ".012") as ifile:
    for line in ifile:
        genoList.append([int(x) for x in line.strip("\n").split("\t")][1:])
geno = np.array(genoList)
del(genoList)
print("Read " + str(geno.shape[1]) + " genotypes for " + str(geno.shape[0]) + " samples")
print("------------------------------------------------")

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
    print("Using prior sex assignments:")
    print("+----------------------------------------------+")
    sexes = pd.read_csv(args.sexes, delimiter="\t", header=None, names=("ind", "sex"))
    sexes = dict(zip(sexes["ind"], sexes["sex"]))
    sexAssign = [0 if sexes[x] == "homogametic" else 1 for x in ind]
    for s in sexes:
        print(str("| " + s).ljust(28, ".") + str(sexes[s] + " |").rjust(20, "."))
    print("+----------------------------------------------+")

else:
    sexes, sexAssign = None, None
    print("No prior sex assignment.")

#########################################
# A function to distinguish uni- vs bi-modal distributions using K-means clustering and a likelihood ratio test
def k2_lrt(dat, init_="k-means++", alpha=float(args.alpha)):
    # init_ is a vector of 0 and 1, describing the a priori groups for initialisation
    # Perform K=2 clustering:
    dat = np.array(dat)
    dat = ma.array(dat, mask=[np.isnan(dat)])

    # Keep indices before removing masked values:
    idx = [x for x in range(len(dat))]
    maidx = ma.array(idx, mask=ma.getmask(dat))
    maidx = maidx[~maidx.mask]  # we keep the idx of unmasked values

    # Remove masked values:
    dat = np.array(dat[~dat.mask])
    dat = dat.reshape(-1, 1)

    # Prepare the assignment list:
    assign = np.full([len(idx)], np.nan)

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
    pval = chi2.sf(LRT, 1 + len(dat))

    if len(dat_0) == 1 or len(dat_1) == 1:  # A one-sample cluster makes no sense (it has no variance)
        pval = 1

    if pval <= alpha:
        isBimodal = 1
        # We polarise the clusters: 0 for low het, 1 for high het
        if mu_2_0 >= mu_2_1:
            assignList = [1 if x == 0 else 0 for x in kmeans.labels_]
        else:
            assignList = [x for x in kmeans.labels_]

    else:
        isBimodal = 0
        assignList = [0] * len(maidx)

    for i, j in enumerate(maidx):
        assign[j] = assignList[i]

    return ([assign, pval, isBimodal])

#########################################
# Compute heterozygosity on blocks of SNPs:
def windowHet_v0(geno, posdf, scaf, wSize=50, wStep=None):

    if not wStep:
        wStep = wSize

    these_idx = [x for x in posdf.loc[posdf["chrom"] == scaf].index]
    these_pos = np.array([x for x in posdf["pos"].loc[these_idx]])
    these_geno = geno[:, these_idx]

    # Bases with "0" and "2" are equally monomorphic. Here we set them all to "0":
    these_geno[these_geno == 2] = 0

    nGeno = these_geno.shape[1]
    print(nGeno)
    nWindows = math.ceil((nGeno - wSize) / wStep) + 1

    these_means = np.full([nInd, nWindows], np.nan)
    wCenterPos = []
    start = 0
    stop = wSize
    w = 0

    while stop <= (nGeno + wStep):
        wData = these_geno[:, start:stop]

        # per-SNP reference allele frequency:
        masked = ma.array(wData, mask=[wData == -1])

        nChroms, nRef = [], []
        for x in range(masked.shape[1]):
            nChroms.append(2 * masked[:, x].count(axis=0))
            nRef.append(np.where(masked[:, x] == 0)[0].shape[0] * 2 + np.where(masked[:, x] == 1)[0].shape[0])
        raf = np.array([2 * x * (1 - x) for x in [nRef[i] / c for i, c in enumerate(nChroms)]])

        # Fold homozygotes:
        masked[masked == 2] = 0

        # Compute adjusted heterozygosity taking into account missing sites:
        F = []
        for h in range(masked.shape[0]):
            this_het = masked[h, :]
            this_raf = ma.array(raf, mask=ma.getmask(this_het))
            F.append(1 - (this_het.mean() / this_raf.mean()))
        these_means[:, w] = ma.array(F, mask=[F is np.nan])

        if stop >= nGeno:
            realStop = nGeno - 1
        else:
            realStop = stop
        wCenterPos.append(np.mean([these_pos[start], these_pos[realStop]]))

        start += wStep
        stop += wStep
        w += 1
    return ([these_means, wCenterPos])


def windowHet(geno, posdf, scaf, wSize=50, wStep=None):
    if not wStep:
        wStep = wSize

    these_idx = [x for x in posdf.loc[posdf["chrom"] == scaf].index]
    these_pos = np.array([x for x in posdf["pos"].loc[these_idx]])
    these_geno = geno[:, these_idx]

    nGeno = these_geno.shape[1]

    nWindows = math.ceil((nGeno - wSize) / wStep) + 1

    these_means = ma.empty([nInd, nWindows])
    wCenter = []
    wCenterPos = []
    start = 0
    stop = wSize
    w = 0

    #while start <= (nGeno - wSize):
    while stop < (nGeno + wStep):
        wData = these_geno[:, start:stop]

        # per-SNP reference allele frequency:
        masked = ma.array(wData, mask=[wData == -1])

        nChroms, nRef = [], []
        for x in range(masked.shape[1]):
            nChroms.append(2 * masked[:, x].count(axis=0))
            nRef.append(np.where(masked[:, x] == 0)[0].shape[0] * 2 + np.where(masked[:, x] == 1)[0].shape[0])

        raf = np.array([2 * x * (1 - x) for x in [nRef[i] / c for i, c in enumerate(nChroms)]])

        # Fold homozygotes:
        masked[masked == 2] = 0
        observed_Het = masked.mean(axis=1)

        F = ma.empty(masked.shape[0])
        for h in range(masked.shape[0]):
            this_het = masked[h, :]
            this_raf = ma.array(raf, mask=ma.getmask(this_het))
            F[h] = 1 - (this_het.mean() / this_raf.mean())

        F = ma.array(F, mask=[np.isnan(x) for x in F])
        these_means[:, w] = F

        if stop >= nGeno:
            realStop = nGeno - 1
        else:
            realStop = stop
        wCenterPos.append(np.mean([these_pos[start], these_pos[realStop]]))

        start += wStep
        stop += wStep
        w += 1

    centers = [x for x in np.array(wCenterPos)[~np.all(these_means.mask, axis=0)]]
    these_means = these_means[:, ~np.all(these_means.mask, axis=0)]

    return ([these_means, centers])

#########################################
# Classify blocks in uni- or bimodal:
def classifyWindows(meanArray, init=None, alpha=float(args.alpha)):

    assignArray = np.full([meanArray.shape[0], meanArray.shape[1]], np.nan)
    pvals, whichBimodals = [], []
    for x in range(meanArray.shape[1]):
        assign, pval, isBimodal = k2_lrt(meanArray[:, x], init_=init, alpha=alpha)
        assignArray[:, x] = assign
        pvals.append(pval)
        whichBimodals.append(isBimodal)

    return ([assignArray, pvals, whichBimodals])

#########################################
# Produce a scaffold-level plot: ## NOTE: It would be nice to produce a genome-level summary plot ####
def makePlot(meanArray, centers, isBimodal, scaf, gff):

    plt.rcParams['figure.figsize'] = [18, 6]

    if gff is not None and scaf is not None:
        # Draw mRNA bands:
        this_gff = gff.loc[gff["scaffold"] == scaf]
        for i in range(this_gff.shape[0]):
            plt.axvspan(this_gff.iloc[i]["start"], this_gff.iloc[i]["end"], color="k", alpha=.1)

    for x in range(meanArray.shape[1]):
        if isBimodal == 1:
            col = ["#61E294" if a == 1 else "#FFC857" for a in assign]
        else:
            col = "#62929E"
        plt.scatter([centers[x]] * meanArray.shape[0], meanArray[:, x], c=col, s=2)

    plt.show()

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
    pval = chi2.sf(LRT, 3 * meanArray.shape[1] - 1 + len(groups)*meanArray.shape[1])
    # One degree of freedom for each additional mean and sd, and one per a priori individual assignment

    return (pval)

########################################################################################################################
# PROCEDURE 1: we process only one scaffold

if args.chrom:

    meanArray, centers = windowHet(geno, posdf, scaf=args.chrom, wSize=wSize, wStep=wStep)
    #print(meanArray)
    if not args.sexes:
        init = "k-means++"
    else:
        init = sexAssign
    assignArray, pvals, whichBimodals = classifyWindows(meanArray=meanArray, init=init)

    if not args.sexes:
        groups = [int(round(x)) for x in np.mean(assignArray, axis=1)]
    else:
        groups = sexAssign
    scaf_pval = ScafLRT(groups, meanArray)

    print(scaf_pval)

########################################################################################################################
# PROCEDURE 2: we process all scaffolds
else:

    scafIsBimodal, groups, meanArrays, assignArrays = [], [], [], []

    if not args.sexes:
        init = "k-means++"
    else:
        init = sexAssign

    # First pass: we determine grouping per scaffold:
    for scaf in scaffolds:

        print(str("Processing " + scaf).ljust(25, " ") + "| Computing heterozygosity |", end="")
        meanArray, centers = windowHet(geno, posdf, scaf=scaf, wSize=wSize, wStep=wStep)
        #print(str("Processing " + scaf).ljust(28, " ") + "| " + str("Classifying SNP blocks").ljust(35, " ") + "|")
        print(" Classifying SNP blocks |", end="")
        assignArray, pvals, whichBimodals = classifyWindows(meanArray=meanArray, init=init)
        these_groups = [x for x in np.mean(assignArray, axis=1)]
        print(str(" " + str(sum(whichBimodals)) + " blocks are bimodal out of " + str(len(whichBimodals))).ljust(35, " ") + "|")
        groups.append(these_groups)
        meanArrays.append(meanArray)
        assignArrays.append(assignArray)

        # If any block is bimodal, we retain the scaffold as potentially testable:
        if not all(x == 0 for x in np.mean(assignArray, axis=1)):
            scafIsBimodal.append(True)
        else:
            scafIsBimodal.append(False)

    #for i, x in enumerate(scafIsBimodal):
    #    if x:
    #        print(assignArrays[i])

    all_assignments = np.concatenate([assignArrays[i] for i, x in enumerate(scafIsBimodal) if x], axis=1)
    allAssignList = []
    for a in range(all_assignments.shape[1]):
        allAssignList.append([x for x in all_assignments[:,a]])
    #c = Counter(str(elem) for elem in allAssignList)
    #for x in c:
    #    print(c[x], x)

    # Initialize list
    countDict = {}

    # Use Append through Iteration
    for elem in allAssignList:
        countDict.setdefault(tuple(elem), list()).append(1)
    for k, v in countDict.items():
        countDict[k] = sum(v)

    # Print Result
    sortedCounts = sorted(zip([x for x in countDict], [countDict[x] for x in countDict]), key=lambda x: x[1], reverse=True)
    sortedVariableCounts = []
    for c in sortedCounts:
        if not all(x == 0 for x in [y for y in c[0] if y == y]):
            #print(str(c[1]) + " :\t" + str(c[0]))
            sortedVariableCounts.append(c)

    majorityGrouping = [x for x in sortedVariableCounts[0][0]]
    asWords = ["heterogametic" if x == 0 else "homogametic" for x in majorityGrouping]
    consensusAssign = dict(zip(ind, asWords))

    #########################################
    # Load the sex assignment file:
    if args.sexes:
        print("Prior sex assignment".ljust(47, " ") + "| Consensus")
        print("+------------------------------------------------------------------+")
        for s in consensusAssign:
            print(str("| " + s).ljust(28, ".") + str(sexes[s] + " |").rjust(20, ".") + (" " + consensusAssign[s]).ljust(19, " ") + "|")
        print("+------------------------------------------------------------------+")

    else:
        print("Consensus sex assignments")
        print("+----------------------------------------------+")
        for s in consensusAssign:
            print(str("| " + s).ljust(28, ".") + str(consensusAssign[s] + " |").rjust(20, "."))
        print("+----------------------------------------------+")

    with open(args.out + ".sex", "w") as ofile:
        for s in consensusAssign:
            ofile.write(s + "\t" + consensusAssign[s] + "\n")

    #########################################
    # Second pass: we evaluate the grouping

    if args.sexes:
        scaf_pvals_prior, scaf_pvals_denovo = [], []
        prior_groups = sexAssign
        denovo_groups = majorityGrouping

        for i, scaf in enumerate(scaffolds):
            print(str("Processing " + scaf).ljust(25, " ") + "| Computing scaffold-level p-value |", end="")
            scaf_pval_prior = ScafLRT(prior_groups, meanArrays[i])
            scaf_pvals_prior.append(scaf_pval_prior)

            scaf_pval_denovo = ScafLRT(majorityGrouping, meanArrays[i])
            scaf_pvals_denovo.append(scaf_pval_denovo)

            if scaf_pval_prior <= args.alpha:
                print(str(" A priori group p-val = " + str(round(scaf_pval_prior, 4)) + "*").ljust(30, " ") + "|", end = "")
            else:
                print(str(" A priori group p-val = " + str(round(scaf_pval_prior, 4))).ljust(30, " ") + "|", end = "")

            if scaf_pval_denovo <= args.alpha:
                print(str(" De novo group p-val = " + str(round(scaf_pval_denovo, 4)) + "*").ljust(30, " ") + "|", end = "\n")
            else:
                print(str(" De novo group p-val = " + str(round(scaf_pval_denovo, 4))).ljust(30, " ") + "|", end = "\n")
