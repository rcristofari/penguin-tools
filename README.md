# penguin-tools

A set of ad-hoc scripts and jupyter notebooks to analyse penguin genomic data. For each tool, see `--help` for the detail of arguments.

`gff2fasta.py` extracts all coding sequences from a reference fasta file, based on gff annotations. It will concatenate the CDS together, adjusting the phasing and strand direction as needed. It can also output the amino-acid sequence with `--type AA`.

`polariseVCF.py` replaces the REF and ALT bases in a VCF file by the ancestral and derived bases, taken from an ancestral reconstruction of the same genome (which needs to share the same coordinates - typically a `baseml`reconstruction. All the matching infromation (genotypes, likelihoods...) follow.

`updateFasta.py` updates a reference genome to include the variants discovered in a specific population. It can replace the reference base with the IUPAC ambiguity code, or the reference / alternative / major / minor allele. If an ancestral reconstruction is provided, ambiguities in the reference genome will be set to the ancestral base, so that only the polymorphism actually present in the variant population will be represented in the final sequence. If the ancestral base is itself ambiguous or undefined, ambiguities are masked.

`markovMask.py` identifies regions of a genome that appear misassembled based on sequencing depth of a panel of individuals. The rationale is that summed sequencing depth for a panel of individuals is expected to follow a negative binomial distribution, and that regions that deviate from this expectation are likely to be collapsed paralogs (over-sequenced) or highly repetitive, hard-to-align DNA (under-sequenced). We fit a negative binomial distribution to the center of the distribution based on quantile-quantile least square regression, and two genetic "flanking" negative binomials to account for high and low coverage areas. We then use a hidden Markov model to identify state changes in the assembly. Finally, outside of outlier blocks, we remove individual sites that seem under- or over-sequenced. The result is a BED file with regions identified as problematic, and a sites file listing files that pass filters.
