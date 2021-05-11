import argparse

parser = argparse.ArgumentParser(description="Identify regions of abnormally high or low sequencing depth, based on outlier regions in an expected negative binomial distribution of sequencing depth. This is based on a hidden Markov model allowing for normal, high or low coverage areas")
parser.add_argument('--idepth', help='Individual mean depth (from vcftools --depth)')
parser.add_argument('--gdepth', help='Per-genotype depth (from vcftools --geno-depth)')
parser.add_argument('--normdepth', help='Pre-normalised depths. If this script has already been run once, this can be provided instead of --idepth and --gdepth')
parser.add_argument('--out', help='Basename of the output BED file', default="depth_mask")
parser.add_argument('--separate', action='store_true',  help='Keep high and low coverage areas in two separate BED file')
parser.add_argument('--quantiles', help='Lower and upper quantiles used for the fit', default=".25,.75")
parser.add_argument('--thin', help='Down-sampling factor for the negative binomial fit (1 for all sites)', default=1000)
parser.add_argument('--trans', help="State transition probability", default=1e-4)
parser.add_argument('--siteThreshold', help="After region masking, remove individual sites below this probability level", default=.01)
parser.add_argument('--chrom', help="Restrict analysis to this scaffold (comma-separated list)", default=None)

args = parser.parse_args()

args.quantiles = [float(x.strip()) for x in args.quantiles.split(",")]
args.trans = float(args.trans)

# Parse the chromosome list argument:
chromlist = None
if args.chrom:
    chromarg = args.chrom.split(",")
    if len(chromarg) > 1:
        chromlist = [x.strip() for x in chromarg]
        args.chrom = None

print("--------------------------------------------------------------------------+")
if args.normdepth:
    print("Normalised depth input file: " + args.normdepth)
else:
    print("Individual depths file: " + args.idepth)
    print("Genotype depths file: " + args.gdepth)
if args.separate:
    print("Low coverage area output: " + args.out + ".lowcov.bed")
    print("High coverage area output: " + args.out + ".highcov.bed")
else:
    print("Outlier coverage area output: " + args.out + ".bed")
print("Outside outlier blocks, filter outliers at a probability level of " + str(args.siteThreshold))
print("Kept sites output: " + args.out + ".pass.sites")
print("--------------------------------------------------------------------------+")
print("Fit parameters:")
print("Retaining %s %% of all sites for negative binomial fit" % str(100/args.thin))
print("Fitting depths between quantiles %s%% and %s%%" % (str(args.quantiles[0]*100), str(args.quantiles[1]*100)))
print("HMM symmetric state transition probability: " + str(args.trans))
print("--------------------------------------------------------------------------+")

import pandas as pd
import numpy as np
from scipy.stats import nbinom
import linecache, random, time

if not args.normdepth:
    # Load in the individual mean depth:
    idepth = pd.read_csv(args.idepth, delimiter="\t")
    idepthDict = dict(zip(idepth["INDV"], [float(x) for x in idepth["MEAN_DEPTH"]]))
    meanDepth = np.mean([float(x) for x in idepth["MEAN_DEPTH"]])

    # Normalise depth per individual mean depth:
    with open(args.gdepth) as ifile:
        with open(args.gdepth.strip(".gdepth") + ".normdepth", "w") as ofile:
            names = next(ifile).strip("\n").split("\t")[2:]
            for line in ifile:
                line = line.strip("\n").split("\t")
                chrom = str(line[0])
                pos = str(line[1])
                sDepth = str(int(round(np.sum([(float(l) / idepthDict[names[i]]) for i, l in enumerate(line[2:])])*meanDepth)))
                ofile.write("\t".join([chrom, pos, sDepth]) + "\n")

    args.normdepth = args.gdepth.strip(".gdepth") + ".normdepth"


########################################################################################################################
# A function to fit the negative binomial quantiles to the observed quantiles:
def gridFit(obs, p=[0, 1], n=[0, 500], domain=[.25, .75], thres=.1):
    q = np.linspace(domain[0], domain[1], 100)

    def ssq(obs, n, p):
        exp = nbinom.ppf(q, n, p)
        ssq = np.sum([(x - exp[i]) ** 2 for i, x in enumerate(obs)])
        return (ssq)

    previous, this = (np.mean(n), np.mean(p)), (0.001, 0.001)
    N = np.linspace(n[0], n[1], 100)
    P = np.linspace(p[0], p[1], 100)

    iter = 0
    while abs(ssq(obs, previous[0], previous[1]) - ssq(obs, this[0], this[1])) > thres:
        iter += 1
        print(str("Iteration # %s" % str(iter)).ljust(15, " ") + "|", end="")
        previous = this[:]
        dist = np.full((100, 100), np.nan)
        for i, ni in enumerate(N):
            for j, pj in enumerate(P):
                d = ssq(obs, ni, pj)
                dist[i, j] = d

        nId = np.where(dist == np.nanmin(dist))[0]
        pId = np.where(dist == np.nanmin(dist))[1]
        this = (np.mean(N[nId]), np.mean(P[pId]))
        nMin, nMax = N[[nId[0] - 10 if nId[0] - 10 > 0 else 0]][0], N[[nId[-1] + 10 if nId[-1] + 10 < 100 else 99]][0]
        pMin, pMax = P[[pId[0] - 10 if pId[0] - 10 > 0 else 0]][0], P[[pId[-1] + 10 if pId[-1] + 10 < 100 else 99]][0]
        # Adjust the edges:
        if pMin == min(P):
            pMin = min(P) * 0.8
        if pMin == 0:
            pMin = 0.0001
        if pMax == max(P):
            pMax = max(P) * 1.5
        if pMax > 1:
            pMax = 1
        if nMin == min(N):
            nMin = min(N) * 0.8
        if nMax == max(N):
            nMax = max(N) * 1.2
        N = np.linspace(nMin, nMax, 100)
        P = np.linspace(pMin, pMax, 100)
        print(str(" Current parameter state %s" % str(round(this[0], 3))).ljust(34, " ") + ";", end="")
        print(str(" %s" % str(round(this[1], 3))).ljust(8, " ") + "|", end="")
        print(str(" SSE = %s" % str(round(np.nanmin(dist), 3))).ljust(14, " ") + "|", end="\n")

    return (this)

########################################################################################################################
# A function to fit the negative binomials by two moments:
def gridOptimFlanks(fitparams, quantiles=(.05, .5), toquantiles=(.5, .99), thres=.1, p=[0, 1], n=[0, 500]):
    # The central distribution:
    fit = nbinom(fitparams[0], fitparams[1])

    # The squared difference between the quantiles:
    def ssq(n, p, quantiles, toquantiles):
        this_fit = nbinom(n, p)
        return (np.sum([(this_fit.ppf(quantiles[i]) - fit.ppf(toquantiles[i])) ** 2 for i in range(len(quantiles))]))

    previous, this = (np.mean(n), np.mean(p)), (0.001, 0.001)
    N = np.linspace(n[0], n[1], 100)
    P = np.linspace(p[0], p[1], 100)

    iter = 0
    while abs(ssq(previous[0], previous[1], quantiles, toquantiles) - ssq(this[0], this[1], quantiles,
                                                                          toquantiles)) > thres:
        iter += 1
        print(str("Iteration # %s" % str(iter)).ljust(15, " ") + "|", end="")
        previous = this[:]
        dist = np.full((100, 100), np.nan)
        for i, ni in enumerate(N):
            for j, pj in enumerate(P):
                d = ssq(ni, pj, quantiles, toquantiles)
                dist[i, j] = d

        np.where(dist == np.nanmin(dist))
        nId = np.where(dist == np.nanmin(dist))[0]
        pId = np.where(dist == np.nanmin(dist))[1]
        this = (np.mean(N[nId]), np.mean(P[pId]))
        nMin, nMax = N[[nId[0] - 10 if nId[0] - 10 > 0 else 0]][0], N[[nId[-1] + 10 if nId[-1] + 10 < 100 else 99]][0]
        pMin, pMax = P[[pId[0] - 10 if pId[0] - 10 > 0 else 0]][0], P[[pId[-1] + 10 if pId[-1] + 10 < 100 else 99]][0]
        # Adjust the edges:
        if pMin == min(P):
            pMin = min(P) * 0.8
        if pMin == 0:
            pMin = 0.0001
        if pMax == max(P):
            pMax = max(P) * 1.5
        if pMax > 1:
            pMax = 1
        if nMin == min(N):
            nMin = min(N) * 0.8
        if nMax == max(N):
            nMax = max(N) * 1.2
        N = np.linspace(nMin, nMax, 100)
        P = np.linspace(pMin, pMax, 100)
        print(str(" Current parameter state %s" % str(round(this[0], 3))).ljust(34, " ") + ";", end="")
        print(str(" %s" % str(round(this[1], 3))).ljust(8, " ") + "|", end="")
        print(str(" SSE = %s" % str(round(np.nanmin(dist), 3))).ljust(14, " ") + "|", end="\n")

    return (this)

########################################################################################################################
# A function to find the best Hidden Markov path using the Viterbi algorithm:
def viterbi(y, A, B, Pi=None, scaf=None):
    # Cardinality of the state space
    start_time = time.time()
    K = A.shape[0]
    # Initialize the priors with default (uniform dist) if not given by caller
    Pi = Pi if Pi is not None else np.full(K, 1 / K)
    T = len(y)
    T1 = np.empty((K, T), 'd')
    T2 = np.empty((K, T), 'B')

    # Initialize the tracking tables from first observation
    T1[:, 0] = np.log(Pi * [b.pmf(y[0]) for b in B])
    T2[:, 0] = 0

    # Iterate through the observations updating the tracking tables
    for i in range(1, T):
        T1[:, i] = np.max(
            T1[:, i - 1] + np.log(A.T * np.array([[b.pmf(y[i]) if b.pmf(y[i]) > 1e-100 else 1e-100 for b in B]]).T), 1)
        # print(y[i])
        # print(np.array([[b.pmf(y[i]) for b in B]]))
        # T1[:, i] = T1[:, i] - np.nanmean(T1[:, i])
        # print(T1[:,i])
        T2[:, i] = np.argmax(T1[:, i - 1] + np.log(A.T), 1)

        if scaf:
            if i % 100000 == 0:
                print(str(scaf + ": processed %s of %s sites" % (i, T)).ljust(65, ".") + str("%.2f s" % (time.time() - start_time)).rjust(10, "."))

    # Build the output, optimal model trajectory
    x = np.empty(T, 'B')
    x[-1] = np.argmax(T1[:, T - 1])
    for i in reversed(range(1, T)):
        x[i - 1] = T2[x[i], i]
    return (x, T1, T2)

########################################################################################################################
# A function to merge together continuous blocks of SNPs in the HMM classification
def mergeBlocks(v, this_scaf, pos):
    start, stop, level = [pos[0]], [], []
    for i, x in enumerate(v[0]):
        if i > 0:
            if x != v[0][i-1]:
                stop.append(pos[i-1])
                level.append(v[0][i-1])
                start.append(pos[i])
    stop.append(pos[-1])
    level.append(v[0][-1])
    return([(this_scaf, s, stop[i], level[i]) for i, s in enumerate(start) if level[i] != 1])

########################################################################################################################
# Implement the workflow:

# Downsample the sites to the specified factor:
nLines = sum(1 for line in open(args.normdepth))
rldepth = [int(linecache.getline(args.normdepth, l).strip("\n").split("\t")[2]) for l in sorted(random.sample([x for x in range(nLines)], round(nLines/float(args.thin))))]

# Fit eh negative binomial distribution
QL, QH = args.quantiles[0], args.quantiles[1]
Qobs = np.quantile(np.array(rldepth), np.linspace(QL, QH, 100))
print("--------------------------------------------------------------------------+")
print("Fitting the central quantiles of the distribution:")
print("Fitting depths between quantiles %s%% and %s%%" % (str(QL*100), str(QH*100)))
print("--------------------------------------------------------------------------+")
start_time = time.time()
n_param, p_param = gridFit(Qobs, domain=[QL, QH], thres=.01)
print("--------------------------------------------------------------------------+")
print("--- Finished in %.2f seconds" % (time.time() - start_time), end="\n\n")
print("Fitting the upper flanking region:")
print("--------------------------------------------------------------------------+")
start_time = time.time()
upperfit = gridOptimFlanks(quantiles=(.1, .5), toquantiles=(.9, .999), fitparams=(n_param, p_param))
print("--------------------------------------------------------------------------+")
print("--- Finished in %.2f seconds" % (time.time() - start_time), end="\n\n")
print("Fitting the lower flanking region:")
print("--------------------------------------------------------------------------+")
start_time = time.time()
lowerfit = gridOptimFlanks(quantiles=(.5, .9), toquantiles=(.001, .1), fitparams=(n_param, p_param))
print("--------------------------------------------------------------------------+")
print("--- Finished in %.2f seconds" % (time.time() - start_time), end="\n\n")

print("+-------------------------------------------------------------------------+")
print("| READY TO SCAN SCAFFOLDS                                                 |")
print("+-------------------------------------------------------------------------+")
# Load chromosome information:
chrom = []
with open(args.normdepth) as ifile:
    for line in ifile:
        chrom.append(line.split("\t")[0])
chromdf = pd.DataFrame(chrom)
del(chrom)

if not chromlist:
    chromlist = sorted([x for x in set(chromdf.iloc[:,0])])

fitHigh = nbinom(upperfit[0], upperfit[1])
fit = nbinom(n_param, p_param)
fitLow = nbinom(lowerfit[0], lowerfit[1])

A = np.array([[1-args.trans, args.trans, args.trans], [args.trans, 1-args.trans, args.trans], [args.trans, args.trans, 1-args.trans]])
B = np.array([nbinom(lowerfit[0], lowerfit[1]), nbinom(n_param, p_param), nbinom(upperfit[0], upperfit[1])])

if args.chrom:
    print("Fitting Hidden Markov Model to " + args.chrom)
    depth = [int(linecache.getline(args.normdepth, i).strip("\n").split("\t")[2]) for i in [x + 1 for x in chromdf.loc[chromdf[0] == args.chrom].index]]
    pos = [int(linecache.getline(args.normdepth, i).strip("\n").split("\t")[1]) for i in [x + 1 for x in chromdf.loc[chromdf[0] == args.chrom].index]]

    hmm = viterbi([int(round(x)) for x in depth], A, B, Pi=None, scaf=args.chrom)
    propLost = round(len(np.where(hmm[0] == 1)[0]) * 100 / len(hmm[0]), 1)
    propLow = round(len(np.where(hmm[0] == 0)[0]) * 100 / len(hmm[0]), 1)
    propHigh = round(len(np.where(hmm[0] == 2)[0]) * 100 / len(hmm[0]), 1)
    print(args.chrom + ": HMM identified %s %% of low-depth and %s %% of high-depth regions" % (propLow, propHigh))
    print(args.chrom + ": after filtering, we retain %s %% of all sites." % propLost)
    dC, dL, dH = np.mean(np.array(depth)[np.where(hmm[0] == 1)[0]]), np.mean(
        np.array(depth)[np.where(hmm[0] == 0)[0]]), np.mean(np.array(depth)[np.where(hmm[0] == 2)[0]])
    print(args.chrom + ": average depth for retained areas: %.1f. " % (dC))
    print(args.chrom + ": average depth for discarded low: %.1f and high: %.1f" % (dL, dH))

    print("Writing out to BED file...", end="")
    bed = mergeBlocks(hmm, args.chrom, pos)
    if not args.separate:
        with open(args.out + "." + args.chrom + ".bed", "w") as ofile:
            for b in bed:
                ofile.write('\t'.join([str(x) for x in b[0:3]]) + '\n')
    else:
        with open(args.out + "." + args.chrom + ".highcov.bed", "w") as highfile, open(args.out + "." + args.chrom + ".lowcov.bed", "w") as lowfile:
            for b in bed:
                if b[3] == 0:
                    lowfile.write('\t'.join([str(x) for x in b[0:3]]) + '\n')
                elif b[3] == 2:
                    highfile.write('\t'.join([str(x) for x in b[0:3]]) + '\n')
    print("done.")

    print("Writing out to kept sites...", end="")
    lowerbound = fit.ppf(float(args.siteThreshold))
    upperbound = fit.ppf(1 - float(args.siteThreshold))
    keptSites = 0
    with open(args.out + "." + args.chrom + ".pass.sites", "a") as sitefile:
        for i, x in enumerate(hmm[0]):
            if x == 1:
                if depth[i] >= lowerbound and depth[i] <= upperbound:
                    sitefile.write(args.chrom + "\t" + str(pos[i]) + "\n")
                    keptSites += 1
    print("done.")
    print("After outlier removal, we retained %.1f %% af all sites in the scaffold" % (keptSites * 100 / len(depth)))

else:
    if not args.separate:
        file = open(args.out + ".bed", "w")
        file.close()
    else:
        file = open(args.out + ".highcov.bed", "w")
        file.close()
        file = open(args.out + ".lowcov.bed", "w")
        file.close()

    # Initialise the kept positions file:
    file = open(args.out + ".pass.sites", "w")
    file.close()

    # Store times taken for 10,000 sites blocks
    sitesDone = 0
    totalSites = len([x for x in chromdf[0] if x in chromlist])
    timePerSite = []

    for idc, c in enumerate(chromlist):
        scaf_time = time.time()
        print(str("Fitting Hidden Markov Model to " + c).ljust(50, " ") + str("(%s / %s)" % (idc+1, len(chromlist))).rjust(25, " "))
        depth = [int(linecache.getline(args.normdepth, i).strip("\n").split("\t")[2]) for i in [x + 1 for x in chromdf.loc[chromdf[0] == c].index]]
        pos = [int(linecache.getline(args.normdepth, i).strip("\n").split("\t")[1]) for i in[x + 1 for x in chromdf.loc[chromdf[0] == c].index]]
        hmm = viterbi([int(round(x)) for x in depth], A, B, Pi=None, scaf=c)
        propLost = round(len(np.where(hmm[0] == 1)[0])*100 / len(hmm[0]), 1)
        propLow = round(len(np.where(hmm[0] == 0)[0])*100 / len(hmm[0]), 1)
        propHigh = round(len(np.where(hmm[0] == 2)[0])*100 / len(hmm[0]), 1)
        print(c + ": HMM identified %s %% of low-depth and %s %% of high-depth regions" % (propLow, propHigh))
        print(c + ": after block filtering, we retain %s %% of all sites." % propLost)
        dC, dL, dH = np.mean(np.array(depth)[np.where(hmm[0] == 1)[0]]), np.mean(np.array(depth)[np.where(hmm[0] == 0)[0]]), np.mean(np.array(depth)[np.where(hmm[0] == 2)[0]])
        print(c + ": average depth for retained areas: %.1f. " % (dC))
        print(c + ": average depth for discarded low: %.1f and high: %.1f" % (dL, dH))

        print("Writing out to BED file...", end="")
        bed = mergeBlocks(hmm, c, pos)
        if not args.separate:
            with open(args.out + ".bed", "a") as ofile:
                for b in bed:
                    ofile.write('\t'.join([str(x) for x in b[0:3]]) + '\n')
        else:
            with open(args.out + ".highcov.bed", "a") as highfile, open(args.out + ".lowcov.bed", "a") as lowfile:
                for b in bed:
                    if b[3] == 0:
                        lowfile.write('\t'.join([str(x) for x in b[0:3]]) + '\n')
                    elif b[3] == 2:
                        highfile.write('\t'.join([str(x) for x in b[0:3]]) + '\n')
        print("done.")

        print("Writing out to kept sites...", end="")
        lowerbound = fit.ppf(float(args.siteThreshold))
        upperbound = fit.ppf(1 - float(args.siteThreshold))
        keptSites = 0
        with open(args.out + ".pass.sites", "a") as sitefile:
            for i, x in enumerate(hmm[0]):
                if x == 1:
                    if depth[i] >= lowerbound and depth[i] <= upperbound:
                        sitefile.write(c + "\t" + str(pos[i]) + "\n")
                        keptSites += 1
        print("done.")
        print("After outlier removal, we retained %.1f %% af all sites in the scaffold" % (keptSites*100/len(depth)))

        timePerSite.append((time.time() - scaf_time)/len(depth))
        sitesDone += len(depth)
        remainingTime = (totalSites-sitesDone)*np.mean(timePerSite) / 3600
        print("Estimated time to completion: %.2f hours" % remainingTime)
        print("+-------------------------------------------------------------------------+")

