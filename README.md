# Chloranthus-sessilifolius-genome
### Genome assembly  
  
### Genome annotation  

### WGD analyses 
all analyses were performed by WGDI: https://github.com/SunPengChuan/wgdi
### Orthologs inferring
```
source sonicparanoid/bin/activate
sonicparanoid -i ./input -o ./output -m fast -t 40  #All the protein sequences were list in the input folder
```
### Orthologs alignment

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
### DiscoVista analyses (for both SSCG and OSCG)

### STAG analyses (for LCG)

### ASTRAL-pro analyses (for LCG)
```
java -D"java.library.path=/data/01/user106/software/astral/A-pro-paper/ASTRAL-MP/lib" -jar /data/01/user106/software/astral/A-pro-paper/ASTRAL-MP/astral.1.1.2.jar -i LCGs.tre -a Ind2Sp.list -o LCG.ASTRAL-pro.tre
# The Ind2Sp.list like this as one species can only retain up to 5 sequences for each cluster:
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

