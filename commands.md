# Commands used in Steinworth_et_al_2020-CnidarianHox
1\. Use  hmm2aln.pl to identify Hox/ParaHox and related homeodomains from translated transcriptomes and protein model files of selected cnidarian datasets.

requires: `HMMer` (http://hmmer.org/) and `hmm2aln.pl` (https://github.com/josephryan/hmm2aln.pl)

 `./hmm2aln.pl --hmm=hd60.hmm --name=HD --fasta_dir=02-RENAMED_DATA --threads=40 --nofillcnf=nofill.hox.conf > cnid_hox_plus.fa`

_Combine all sequences in cnid_hox_plus.fa with bilaterian and known cnidarian homeoboxes to create file, all_hox_plus.fa_

2\. Remove sequences with 5 or more gaps using custom script
`./nogaps.py all_hox_plus.fa`

_Add back in sequences cloned from Cassiopea xamachana to file all_hox_plus.fa_wholeSeqs and call this file all_hox_plus.fa_

3\. Generate an initial phylogenetic tree using resulting alignment from `hmm2aln.pl`

requires IQ-tree (http://www.iqtree.org/)

`iqtree-omp -s final_all_hox_plus.fa -nt AUTO -bb 1000 -m LG -pre MLtree_withgaps > iq.out 2> iq.err`

4\. Using resulting tree and alignment, prune non-Hox/ParaHox genes using make subalignment
(requires `make_subalignment_fasta` (https://github.com/josephryan/make_subalignment_fasta)
`./make_subalignment_fasta --tree=MLtree_withgaps.treefile --aln=final_all_hox_plus.fa --root=Anthopleura_elegantissima_45195 --pre=Nvec`

5\. Using the resulting alignment from make_subalignment run the following ML trees 

  a\. RAXML with 25 starting parsimony trees
  `raxmlHPC-PTHREADS-SSE3 -T 25 -p 1234 -# 25 -m PROTGAMMALG -s subalign -n raxML_mp`

  b\. RAXML with 25 random starting trees
  `raxmlHPC-PTHREADS-SSE3 -T 25 -d -p 1234 -# 25 -m PROTGAMMALG -s subalign -n raxML_rt`

  c\. IQTREE
  `iqtree -m LG+G4 -s subalign -pre iqtree -bb 1000`

6\. Evaluate the likelihood scores of the IQTREE using RAxML (for apples-to-apples comparison of likelihoods)
`raxmlHPC-SSE3 -f e -m PROTGAMMALG -t iqtree.treefile -s subalign -n iqtree_raxml`

_Compare Final GAMMA-based score of RAXML trees and IQTREE tree_ 

7\. For the best tree (in our case this was the RAxML maximum parsimony starting tree) run and apply bootstraps

`raxmlHPC -m PROTGAMMALG -s subalign -p 12345 -x 12345 -# 100 -n raxml_mp_best `

`raxmlHPC -m PROTGAMMALG -p 12345 -f b -t RAxML_bestTree.raxml_mp -z RAxML_bootstrap.raxml_mp_best -n T15`

8\. Run Bayesian tree

`fasta2phy.pl subalign > sub.phy
phy2bayesnex.pl sub.phy > hox.nex`

_Paste the following execution block into hox.nex:_
`mcmcp ngen=10000000 samplefreq=10000 mcmcdiagn=yes stoprule=yes stopval=0.01       nruns=2 nchains=5 savebrlens=yes; mcmc; sumt filename=FILE.nex nRuns=2 Relburnin=YES BurninFrac=.25 Contype=Allcompat;).`

`mpirun -np 25 mb hox.nex`

_Run ML trees to calculate final GAMMA-based score for Bayesian trees_

`raxmlHPC-SSE3 -f e -m PROTGAMMALG -t hox.nex.run1.newick -s ../subalign -n bayes_run1_raxml

raxmlHPC-SSE3 -f e -m PROTGAMMALG -t 03_hox.nex.run2.newick -s ../subalign -n bayes_run2_raxml

raxmlHPC-SSE3 -f e -m PROTGAMMALG -t hox.nex.con.newick -s ../subalign -n bayes_con_raxml`

_Compare final GAMMA-based score to RAxML trees and select best tree for main figure`


AU Test
perl make_constraint_trees.pl > iq_script.sh

This script creates constraint trees for every paired combo of the following:
   1. cnidarian + cnidarian
   2. cnidarian + bilat
   3. cnidarian + mixed (mixed = clade like Gbx includes cnid + bilat)

It also creates 4 additional constraint trees:
   1. ax6,ax6a,bilat_hox1
   2. hd065,cdx,xlox
   3. ax1,ax1a,bilat_post
   4. ax1,ax1a,bilat_post,cent

Also prints out iqtree command lines to STDOUT which we directed to iq_script.sh

cat *.treefile > autest.treels

iqtree -s subalign -m LG+G4 -z autest.treels -n 0 -zb 1000 -au

Removing duplicate sequences from final tree
java -jar /usr/local/phyutility/phyutility.jar -pr -in TREENAME -out TREENAME.pruned -names NAME1 NAME2 NAME3

java -jar /usr/local/phyutility/phyutility.jar -pr -in RAxML_bipartitions.T15 -out RAxML_bipartitions.T15.pruned -names Corallium_rubrum_87184, Corallium_rubrum_87185, Corallium_rubrum_17719, Craspedacusta_sowerbyi_16174, Calvadosia_cruxmelitensis_1956, Corallium_rubrum_75060, Haliclystus_sanjuanensis_2068, Craspedacusta_sowerbyi_25502, Craspedacusta_sowerbyi_93900, Corallium_rubrum_89239, Corallium_rubrum_91898, Eunicella_cavolinii_29512, Alatina_alata_18573, Craspedacusta_sowerbyi_43026, Craspedacusta_sowerbyi_85141, Corallium_rubrum_75469, Cxam_transcriptome_73398

java -jar /usr/local/phyutility/phyutility.jar -pr -in RAxML_bipartitions.T15.pruned.renamed -out RAxML_bipartitions.T15.pruned.renamed -names Cxam_transcriptome_40272

^just because I forgot to remove Cxam_transcriptome_40272
