# PLANNED ANALYSES FOR PHYLOGENETIC ANALYSIS OF CNIDARIAN HOX AND HOXL GENES
 Principle Investigator: Joseph Ryan and Bailey Steinworth
 
 Draft or Version Number: v.1.1
 
 Date: 13 September 2022

## LIST OF ABBREVIATIONS

ax1  - Anthox1

ax1a - Anthox1a

ax6  - Anthox6

ax6a - Anthox6a

ax7  - Anthox7

ax8a - Anthox8a

ax8b - Anthox8b

ax9  - Anthox9

*Nvec* - *Nematostella vectensis* (Hexacorallia)

*Aele* - *Anthopleura elegantissima* (Hexacorallia)

*Adig* - *Acropora digitifera* (Hexacorallia)

*Lscu* - *Lobactis scutaria* (Hexacorallia)

*Crub* - *Corallium rubrum* (Octocorallia)

*Ecav* - *Eunicella cavolinii* (Octocorallia)

*Rren* - *Renilla reniformis* (Octocorallia)

*Hsan* - *Haliclystus sanjuanensis* (Staurozoa)

*Ccru* - *Calvadosia cruxmelitensis* (Staurozoa)

*Cfle* - *Chironex fleckeri* (Cubozoa)

*Aala* - *Alatina alata* (Cubozoa)

*Cxam* - *Cassiopea xamachana* (Scyphozoa)

*Cfus* - *Chrysaora fuscescens* (Scyphozoa)

#*Chem* - *Clytia hemisphaerica (Hydrozoa)

*Hvul* - *Hydra vulgaris* (Hydrozoa)

*Csow* - *Craspedacusta sowerbii* (Hydrozoa)

*Hsap* - *Homo sapiens* (Bilateria)

*Bflo* - *Branchiostoma floridae* (Bilateria)

*Dmel* - *Drosophila melanogaster* (Bilateria)

*Tcas* - *Tribolium castaneum* (Bilateria)

*Ctel* - *Capitella teleta* (Bilateria)

*Cgig* - *Crassostrea gigas* (Bilateria)

## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE

### 1.1 _Background Information_

The relationship of cnidarian and bilaterian HOXL class genes is not well understood (Ryan et al, 2006; Ryan et al, 2007; Dubuc et al, 2011).

### 1.2 _Rationale_

Previous studies have been limited to a very few Cnidaria taxa (e.g., Nematostella vectensis, Hydra magnipapillata, Clytia hemisphaerica, and Acropora digitifera) and usually only a few homeoboxes from bilaterian species.

### 1.3 _Objectives_

Accurately classify cnidarian Hox and ParaHox genes from a wide array of cnidarians.

## 2 STUDY DESIGN AND ENDPOINTS

#### 2.1 Build dataset.

We will start with four sets of Hox/ParaHox homeodomains from the curated HomeoDB (Zhong et al, 2008): *Homo sapiens*, *Branchiostoma floridae*, *Drosophila melanogaster*, and *Tribolium castaneum*. We will include the spiralia Hox/ParaHox homeodomains *Capitella teleta* and *Crassostrea gigas* (as classified in Paps et al, 2015 and Zwarycz et al, 2015). We will include the Hox/ParaHox genes from *Nematostella vectensis* and *Acropora digitifera* (as classified in DuBuc et al, 2012). We will include the *Hydra magnipapillata* and *Clytia hemisphaerica* Hox/ParaHox genes (as classified in Chiori et al, 2009). NOTE: HmaCnox4 is listed as pirS39067; we could not find this so are using CAA45911.1. We will use hd60.hmm hidden Markov model from Zwarycz et al 2015 to search transcriptomes from the following:  *Anthopleura elegantissima*, *Lobactis scutaria*, *Corallium rubrum*, *Renilla reniforms*, *Eunicella cavolinii*, *Haliclystus sanjuanensis*, *Calvadosia cruxmelitensis*, *Chironex fleckeri*, *Alatina alata*, *Cassiopea xamachana*, *Chrysaora fuscescens*, and *Craspedacusta sowerbyi*.

This script runs hmmsearch, stockholm2fasta, and some custom code to remove indels and fill end gaps.
```
./hmm2aln.pl --hmm=hd60.hmm --name=HD --fasta_dir=02-RENAMED_DATA --threads=40 --nofillcnf=nofill.hox.conf > cnid_hox_plus.fa
```

We will remove from the resulting alignment any sequences with five or more gaps (8.3%).

#### 2.2 Generate an initial phylogenetic tree

```
iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err
```

#### 2.3 Prune non-Hox/ParaHox genes

This script takes a prefix of a subset of (our ingroup) taxa and will return an alignment of only those genes within the clade descended from the most recent common ancestor of all genes with the specified prefix. We will use Nematostella as the ingroup prefix.
```
./make_subalignment2 --tree=<newick_treefile> --aln=<phylip_alignment> --root=<root_taxa> --pre=<prefix>
```
#### 2.4 Add cloned Cxam sequences and repeat pruning step
After cloning all _C. xamachana_ genes identified and sequencing, we will add back in the cloned sequences for genes that had five or more gaps.
```
raxmlHPC-SSE3.PTHREADS -T 25 -p [random_number] -# 25 -m PROTGAMMA[best-fit_model] -s [alignment_file] -n [name]_mp
```
./make_subalignment2 --tree=<newick_treefile> --aln=<phylip_alignment> --root=<root_taxa> --pre=<prefix>
```
```

#### 2.5  RAXML with 25 starting parsimony trees and 25 random starting trees; the best fit model will be LG (but we will confirm this using one RAXML run with PROTGAMMAAUTO.

```
raxmlHPC-SSE3.PTHREADS -T 25 -p [random_number] -# 25 -m PROTGAMMA[best-fit_model] -s [alignment_file] -n [name]_mp
```
```
raxmlHPC-SSE3.PTHREDAS -T 25 -d -p [random_number] -# 25 -m PROTGAMMA[best-fit_model] -s [alignment_file] -n [name]_rt
```

### 2.6 IQTREE
```
iqtree -m [best-fit_model]+G4 -s [alignment_file] -pre [prefix] -bb 1000
```

#### COMPARE Iqtree and 50 rax trees using rax to report the likelihood values; generate a likelihood score using RAxML for Iq-tree and grep for the likelihood values from RAxML_info files for RAxML runs.

```
raxmlHPC-SSE3 -f e -m PROTGAMMA[best-fit_model] -t [single-gene_tree] -s [alignment_file] -n [output_name]
grep 'Starting final GAMMA-based' *info*
```

### Run bootstraps and apply them to the best tree above
```
raxmlHPC -m PROTGAMMA[best-fit_model] -s [alignment_file] -p 12345 -x 12345 -# 100 -n [name]
```

```
raxmlHPC -m PROTGAMMA[best-fit_model] -p 12345 -f b -t RAxML_bestTree.[name] -z RAxML_bootstrap.[name] -n T15
```

#### 2.7 Mr. Bayes

We will run 5 MrBayes runs with the following command:

```mpirun -np 25 mb hox.nex```

and the following execution block:

```prset aamodelpr = fixed(LG); lset rates = gamma;
mcmcp ngen=10000000 samplefreq=10000 mcmcdiagn=yes stoprule=yes stopval=0.01
      nruns=2 nchains=5 savebrlens=yes;
mcmc;
sumt filename=FILE.nex nRuns=2 Relburnin=YES BurninFrac=.25 Contype=Allcompat;
```

### CHOOSE TREE THAT WILL BE MAIN FIGURE

Since we cannot use Bayesian principles to evaluate our ML trees, we will use ML principles to evaluate the Bayes trees. If Bayes tree has better likelihood score than the ML tree we will report the Bayes tree as our main figure with BS values from above. If our ML tree has a higher likelhood than our Bayes tree we will report the ML tree with Bayesian log likelihood scores.
All trees will be presented supplement.  And differences between Bayes and ML will be discussed.

#### 2.8 AU Test

We will test all sister-group relationships between cnidarian and cnidarian or bilaterian genes, plus the following topologies, with the AU test using IQTREE
```
((ax6a,ax6,Hox1),all,other,clades)
((ax7,ax8a,ax8b,Hox2),all,other,clades)
((cdx,HD065,xlox),all,other,clades)
((ax1,ax1a,bilat-post),all,other,clades)
((ax1,ax1a,bilat-post,bilat-cent),all,other,clades)
```

#### 2.8 How we will report results

* Report clades as-is where best Baysian and best ML tree agree 

* Determine the likelihood score for best ML tree and best Bayesian tree (choose ML tree if Bayesian tree has a polytomy). Report clades from the tree with the highest likelihood, mention the conflict.  Report summaries of both trees as main figures. 


## 3 LITERATURE REFERENCED

Chiori R, Jager M, Denker E, Wincker P, Da Silva C, Le Guyader H, Manuel M, Quéinnec E. Are Hox genes ancestrally involved in axial patterning? Evidence from the hydrozoan Clytia hemisphaerica (Cnidaria). PloS one. 2009 https://doi.org/10.1371/journal.pone.0004231

DuBuc TQ, Ryan JF, Shinzato C, Satoh N, Martindale MQ. Coral comparative genomics reveal expanded Hox cluster in the cnidarian–bilaterian ancestor. Integrative and Comparative Biology. 2012 https://doi.org/10.1093/icb/ics098

Ryan JF, Mazza ME, Pang K, Matus DQ, Baxevanis AD, Martindale MQ, Finnerty JR. Pre-bilaterian origins of the Hox cluster and the Hox code: evidence from the sea anemone, Nematostella vectensis. PloS One. 2007 https://doi.org/10.1371/journal.pone.0000153

Ryan JF, Burton PM, Mazza ME, Kwong GK, Mullikin JC, Finnerty JR. The cnidarian-bilaterian ancestor possessed at least 56 homeoboxes: evidence from the starlet sea anemone, Nematostella vectensis. Genome Biology. 2006 https://doi.org/10.1186/gb-2006-7-7-r64

Paps J, Xu F, Zhang G, Holland PW. Reinforcing the egg-timer: Recruitment of novel Lophotrochozoa homeobox genes to early and late development in the Pacific oyster. Genome Biology and Evolution. 2015 https://doi.org/10.1093/gbe/evv018

Zhong YF, Butts T, Holland PW. HomeoDB: a database of homeobox gene diversity. Evolution & Development. 2008 https://doi.org/10.1111/j.1525-142X.2008.00266.x

Zwarycz AS, Nossa CW, Putnam NH, Ryan JF. Timing and scope of genomic expansion within Annelida: evidence from homeoboxes in the genome of the earthworm Eisenia fetida. Genome Biology and Evolution. 2015 https://doi.org/10.1093/gbe/evv243

## 5 PHYLOTOCOL AMENDMENT HISTORY

### 1.1 (13 September 2022)
Additions:
* Added cloned sequences for Cxam genes that had been removed for having 5 or more gaps (step 2.4)
* Added second pruning step after Cxam genes were added back into alignment (step 2.4)
* For the AU test, instead of only testing a selection of hypotheses, we tested every possible combination of cnidarian clade plus cnidarian or bilaterian clade, in addition to the other select combinations listed (step 2.8).
