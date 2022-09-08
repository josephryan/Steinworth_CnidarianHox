# Commands used in Steinworth-CnidarianHox
1\. Use  hmm2aln.pl to identify Hox/ParaHox and related homeodomains from translated transcriptomes and protein model files of selected cnidarian datasets.

requires: `HMMer` (http://hmmer.org/) and `hmm2aln.pl` (https://github.com/josephryan/hmm2aln.pl)

 `./hmm2aln.pl --hmm=hd60.hmm --name=HD --fasta_dir=transcriptomes --threads=40 --nofillcnf=nofill.hox.conf > cnid_hox_plus.fa`

_Combine all sequences in cnid_hox_plus.fa with bilaterian and known cnidarian homeoboxes to create file, all_hox_plus.fa_

2\. Remove sequences with 5 or more gaps using custom script

 `./nogaps.py all_hox_plus.fa`

_Remove duplicate sequences_

3\. Generate an initial phylogenetic tree using resulting alignment from `hmm2aln.pl`

requires IQ-tree (http://www.iqtree.org/)

 `iqtree-omp -s all_hox_plus.fa_wholeSeqs -nt AUTO -bb 1000 -m LG -pre IQTree_Initial > iq.out 2> iq.err`

4\. Using resulting tree and alignment, prune non-Hox/ParaHox genes using make subalignment
(requires `make_subalignment2` https://github.com/josephryan/Steinworth_CnidarianHox/blob/master/03-SCRIPTS/make_subalignment2)

 `./make_subalignment2 --tree=IQTree_Initial.treefile --aln=all_hox_plus.fa_wholeSeqs --root=C_sowe.0035330 --pre=Nvec > subalign`

_Add sequences cloned from Cassiopea xamachana to the subalignment_

5\. Using the resulting alignment from make_subalignment2 run the following ML trees 

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

 `raxmlHPC -m PROTGAMMALG -s subalign -p 12345 -x 12345 -# 100 -n raxml_mp_best`
 
 `raxmlHPC -m PROTGAMMALG -p 12345 -f b -t RAxML_bestTree.raxml_mp -z RAxML_bootstrap.raxml_mp_best -n T15`

8\. Run Bayesian tree

`fasta2phy.pl subalign > sub.phy`
`phy2bayesnex.pl sub.phy > hox.nex`

_Paste the following execution block into hox.nex:_
`mcmcp ngen=10000000 samplefreq=10000 mcmcdiagn=yes stoprule=yes stopval=0.01       nruns=2 nchains=5 savebrlens=yes; mcmc; sumt filename=FILE.nex nRuns=2 Relburnin=YES BurninFrac=.25 Contype=Allcompat;).`

`mpirun -np 25 mb hox.nex`

_Run ML trees to calculate final GAMMA-based score for Bayesian trees_

`raxmlHPC-SSE3 -f e -m PROTGAMMALG -t hox.nex.run1.newick -s ../subalign -n bayes_run1_raxml`

`raxmlHPC-SSE3 -f e -m PROTGAMMALG -t 03_hox.nex.run2.newick -s ../subalign -n bayes_run2_raxml`

`raxmlHPC-SSE3 -f e -m PROTGAMMALG -t hox.nex.con.newick -s ../subalign -n bayes_con_raxml`

_Compare final GAMMA-based score to RAxML trees and select best tree for main figure_


9\. Run AU Test

 `perl make_constraint_trees.pl > iq_script.sh`

_This script creates constraint trees for every paired combo of the following:_

   _1. cnidarian + cnidarian_
   
   _2. cnidarian + bilat_
   
   _3. cnidarian + mixed (mixed = clade like Gbx includes cnid + bilat)_

_It also creates 4 additional constraint trees:_

   _1. ax6,ax6a,bilat_hox1_
   
   _2. hd065,cdx,xlox_
   
   _3. ax1,ax1a,bilat_post_
   
   _4. ax1,ax1a,bilat_post,cent_

_It also prints out iqtree command lines to STDOUT which we directed to_ `iq_script.sh`

 `cat *.treefile > autest.treels`
 
 `iqtree -s subalign -m LG+G4 -z autest.treels -n 0 -zb 1000 -au`

