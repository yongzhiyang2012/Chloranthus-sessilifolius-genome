# Chloranthus-sessilifolius-genome
### Genome assembly
```
smartdenovo.pl -p prefix -c 1 reads.fa > prefix.mak
wtpre -J 3,000
wtzmo-k 21 -z 10 -Z 16 -U -1 -m 0.1 -A 1,000
wtclp -d 3 -k 300 -m 0.1 -FT
wtlay -w 300 -s 200 -m 0.1 -r 0.95 -c 1
smartdenovo.pl -p prefix -c 1 reads.fa > prefix.mak
make -f prefix.mak
```
### Genome annotation
```
de novo:
genscan Arabidopsis.smat contig1.fa| phase.genscan.output.pl - > contig1.fa.gff3
augustus --species=BUSCO_final.contig.fa.busco-long_xxxx contig1.fa | perl ConvertFormat_augustus.pl - contig1.gff
trainGlimmerHMM contig1.fa  exon.file -d  glimmerhmm_contig1  2>&1 | tee trainGlimmerHMM_contig1.log
glimmerhmm_linux_x86_64 contig1.fa -d glimmerhmm_contig1  -g  -o  tmp/contig1.gff ; perl ConvertFormat_glimmerhmm.pl tmp/contig1.fa glimmerhmm_results/contig1.gff

homolog-based:
GeMoMa GeMoMaPipeline threads=40 t=final_assembly_genome.fa s=own a=reference.gff  g=reference.fa  outdir=XXX  AnnotationFinalizer.r=NO tblastn=false

transcriptome-based:
perl Launch_PASA_pipeline.pl -c final_assembly_genome.fa -C -R -g ./xxx  -T -t Trinity.fasta.clean -u trinity_all.Trinity.fasta  --ALIGNERS blat,gmap --CPU 50

EVM:
perl partition_EVM_inputs.pl --genome ./final_assembly_genome.fa --gene_predictions ./all_abinitio.gff --protein_alignments ./all_homologs.gff --transcript_alignments ./jsl_mydb_pasa.sqlite.pasa_assemblies.gff3 --segmentSize 1000000 --overlapSize 200000 --partition_listing partitions_list.out
perl write_EVM_commands.pl --genome ./final_assembly_genome.fa --weights evm.weights.txt --gene_predictions ./all_abinitio.gff --protein_alignments ./all_homologs.gff --transcript_alignments ./jsl_mydb_pasa.sqlite.pasa_assemblies.gff3 --partitions partitions_list.out --output_file_name evm.out  > commands.list
perl recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file evm.out
perl convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome final_assembly_genome.fa
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3

repeat annotation:
RepeatMasker -e rmblast -pa 60 -nolow -norna -no_is -gff -species Mesangiospermae final.contig.fa

BuildDatabase -name maoli_db final.contig.fa  2>&1 | tee repeatmodeler.log
RepeatModeler -pa 60 -database maoli_db 2>&1 | tee 02.RepeatModeler.log
RepeatMasker -pa 30 -lib ./custom.lib final.contig.fa

gt suffixerator -db final.contig.fa -indexname final.contig -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest  -index final.contig  -similar 90 -vic 10 -seed 20 -seqids yes  -minlenltr 60 -maxlenltr 7000 -mintsd 4 -maxtsd 20 -motifmis 3  > genome.fasta.harvest.scn
ltr_finder Contig1.fa -w 2 > Contig1.fa.out
LTR_retriever -genome final.contig.fa -inharvest genome.fasta.harvest.scn -infinder genome.finder.scn -threads 30
```

### WGD analyses 
all analyses were performed by WGDI: https://github.com/SunPengChuan/wgdi
### Orthologs inferring
```
sonicparanoid:
source sonicparanoid/bin/activate
sonicparanoid -i ./input -o ./output -m fast -t 40  #All the protein sequences were list in the input folder

orthomcl:
diamond makedb --in goodProteins.fasta -d goodProteins.fasta
diamond blastp -d goodProteins.fasta -q goodProteins.fasta -o diamond.blastpout -p 75 --more-sensitive --outfmt 6
orthomclBlastParser diamond.blastpout compliantFasta > similarSequences.txt
perl -p -i -e 's/0\t0/1\t-181/' similarSequences.txt
orthomclInstallSchema orthomcl.config.template
orthomclLoadBlast orthomcl.config.template similarSequences.txt
orthomclPairs orthomcl.config.template orthomcl_pairs.log cleanup=no
orthomclDumpPairsFiles orthomcl.config.template
mcl mclInput --abc -I 1.5 -o mclOutput
mcl mclInput --abc -I 2.0 -o mclOutput_I2.0
orthomclMclToGroups cluster 1 < mclOutput > groups.txt
```
### Orthologs alignment and trim
```
mafft --auto --quiet gene.pep > gene.pep.aln
pal2nal.pl gene.pep.aln gene.cds -output fasta > gene.cds.aln
pxclsq -a -p 0.2 -s gene.cds.aln -o gene.cds.aln-0.2cln
```
### Concatenation or gene tree building (for both SSCG, OSCG, LCG and chloroplast genes)
```
iqtree -s cds.best.fas.gt80 -st DNA -pre cds.best.fas.gt80 -nt 20 -bb 1000 -m MFP -quiet -redo
```
### ASTRAL analyses (for both SSCG and OSCG)
```
java -jar astral.5.7.1.jar -i SSCG.gene.tre -o SSCG.ASTRAL.tre
java -jar astral.5.7.1.jar -i SSCG.gene.tre -q SSCG.ASTRAL.tre -o SSCG.ASTRAL.tre.t8 -t 8
java -jar astral.5.7.1.jar -i OSCG.gene.tre -o OSCG.ASTRAL.tre
java -jar astral.5.7.1.jar -i OSCG.gene.tre -q OSCG.ASTRAL.tre -o OSCG.ASTRAL.tre.t8 -t 8
```
### Densitree analyses (for SSCG)
```
library(ape)
library(Matrix)
library(phybase)
tree_str <- "tree.file"
name <- species.name(tree_str)
name_len <- length(name)
tree_node <- read.tree.nodes(tree_str,name)$nodes
node_len <- nrow(tree_node)
tree_node[,4] <- 1
new_node <- tree.noclock2clock(node_len,tree_node,name_len)
node_height <- node.height(node_len,new_node,name_len)
new_node[,4] <- new_node[,4]*(1/node_height)
final_tree <- write.subtree(node_len,new_node,name,node_len)
cat(final_tree)
cat("
")
```
### DiscoVista analyses (for both SSCG and OSCG)
```
DiscoVista/src/utils/generate_clade-defs.py annotation.txt clade-defs.txt
DiscoVista/src/utils/discoVista.py -c clade-defs.txt -p gene_tree -w newOrders.txt -t 75 -m 1 -o gene_tree/results
```
### STAG analyses (for LCG)
```
stag  SpeciesMap.txt  Gene_Trees
```

### ASTRAL-pro analyses (for LCG)
```
java -D"java.library.path=/data/01/user106/software/astral/A-pro-paper/ASTRAL-MP/lib" -jar /data/01/user106/software/astral/A-pro-paper/ASTRAL-MP/astral.1.1.2.jar -i LCGs.tre -a Ind2Sp.list -o LCG.ASTRAL-pro.tre
# As one species can only retain up to 5 sequences for each cluster, so the Ind2Sp.list looks like this:
# Aco_1 Aco
# Aco_2 Aco
# Aco_3 Aco
# Aco_4 Aco
# Aco_5 Aco
# Ppe_1 Ppe
# Ppe_2 Ppe
# Ppe_3 Ppe
# Ppe_4 Ppe
# Ppe_5 Ppe
# ...
```

### Hybridization analyses (PhyloNetworks running in Julia programming language with SSCG trees)
```
using Distributed
addprocs(50)  #using 50 threads
using PhyloNetworks;
iqtreeCF = readTrees2CF("./SSCG.gene.tre", writeTab=faslse, writeSummary=false)
astraltree = readMultiTopology("./SSCG.ASTRAL.tre")[1]
net0 = snaq!(astraltree,iqtreeCF, hmax=0, filename="net0", seed=1234, runs=100)
net1 = snaq!(net0, iqtreeCF, hmax=1, filename="net1", seed=1234, runs=100)
net2 = snaq!(net1, iqtreeCF, hmax=2, filename="net2", seed=1234, runs=100)
net3 = snaq!(net2, iqtreeCF, hmax=3, filename="net3", seed=1234, runs=100)
```
### ILS analyses (Simulation using Phybase and DendroPy with the SSCG dataset)
#### Phybase (See our previous work: https://github.com/wk8910/ILS_simulation and https://github.com/yongzhiyang2012/Euryale_ferox_and_Ceratophyllum_demersum_genome_analysis)
#### DendroPy (We used the same pipeline that was applied in avian and mammal analyses. https://www.ideals.illinois.edu/bitstream/handle/2142/55627/Simulation%20Procedure.html?sequence=4&isAllowed=y)
```
wget -O generateCoalescentTrees.py https://www.ideals.illinois.edu/bitstream/handle/2142/55627/generateCoalescentTrees.py?sequence=2  
python2 generateCoalescentTrees.py SSCG.ASTRAL.tre 20000 DendroPy.simulated.tre
```
### Transcriptome analyses
```
hisat2-build -p 20 reference.fa reference
hisat2 -x reference -p 20 --dta --summary-file SRRxxxxxxx_log -U SRRxxxxxxx.fastq.gz | samtools view -@ 20 -bh | samtools sort -@ 20 -o SRRxxxxxxx_sorted.bam

stringtie -o SRRxxxxxxx.gtf  -p 20 -A SRRxxxxxxx.tab -G reference.gff3 -B -e SRRxxxxxxx_sorted.bam
stringtie --merge -e -p 10 -G reference.gff3  -o stringtie_merged.gtf mergelist
stringtie -e -B -p 8 -G stringtie_merged.gtf -o SRRxxxxxxx_merged.gtf SRRxxxxxxx_sorted.bam

python /data/01/user119/tools/py/prepDE.py -i sample_list -l 150
```
