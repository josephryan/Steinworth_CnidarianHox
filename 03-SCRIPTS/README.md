# Scripts used in Steinworth et al. 2020

### make_subalignment2
This script takes a prefix of a subset of (our ingroup) taxa and will return an alignment of only those genes within the clade descended from the most recent common ancestor of all genes with the specified prefix. It differs from https://github.com/josephryan/make_subalignment only in that it takes a FASTA file as an input instead of PHYLIP

### nogaps.py
This script removes all sequences with five or more gaps from a FASTA alignment.
