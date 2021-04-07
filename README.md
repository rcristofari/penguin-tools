# penguin-tools

A set of ad-hoc scripts and jupyter notebooks to analyse penguin genomic data

`gff2fasta.py` extracts all coding sequences from a reference fasta file, based on gff annotations. It will concatenate the CDS together, adjusting the phasing and strand direction as needed. It can also output the amino-acid sequence with `--type AA`.

`polariseVCF.py` replaces the REF and ALT bases in a VCF file by the ancestral and derived bases, taken from an ancestral reconstruction of the same genome (which needs to share the same coordinates - typically a `baseml`reconstruction. All the matching infromation (genotypes, likelihoods...) follow.

`updateFasta.py` updates a reference genome to include the variants discovered in a specific population. It can replace the reference base with the IUPAC ambiguity code, or the reference / alternative / major / minor allele. If an ancestral reconstruction is provided, ambiguities in the reference genome will be set to the ancestral base, so that only the polymorphism actually present in the variant population will be represented in the final sequence. If the ancestral base is itself ambiguous or undefined, ambiguities are masked.
