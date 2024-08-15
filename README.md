# Pangenome Graphs - Simulation

This repository contains a section of the original [Pangenome Graphs Workshop](https://github.com/GenomicsAotearoa/Pangenome-Graphs-Workshop) that was removed to reduce the amount of content being covered.

It involves the use of the `wgsim` software to simulate individual-level genomes (with genetic variants) from a single reference genome represented in `.fasta` format. These multiple genomes from individuals of the same species can then be used for testing/benchmarking pangenome graph construction and variant calling methods.

## Genomic data

The data for this exercise is the same as that used in the [Pangenome Graphs Workshop](https://github.com/GenomicsAotearoa/Pangenome-Graphs-Workshop), which can be download via:

```
git clone https://github.com/ZoeYang2020/dataset_for_pg_workshop
```

The data file `5NM.fa` contains the genome sequence for five different assesmblies of *Neisseria meningitidis*, which need to be split into individual genomes. 

We can see the genome names via:

```
grep '^>' 5NM.fa
```

The above command searches the file `5NM.fa` for lines that start ('^') with '>'.

We can split the file into the five separate genomes (using the '>' symbol to identify where each one starts) usign the following commands:

```
awk '/^>/ { gsub(">","",$1)
            FILE=$1 ".fa"
            print ">" $1 >> FILE
            next}
          { print >> FILE }' 5NM.fa   
```

We now have five separate `.fasta` files, one per genome:

```
NC_003112.2.fa
NC_017518.1.fa
NZ_CP007668.1.fa
NZ_CP016880.1.fa
NZ_CP020423.2.fa
```

## Use `wgsim` to simulate individual genomes and variation

Use `wgsim` to simulate 2x150bp NGS data, with an error rate of 0.005. 

```
outputfolder=simulated-genomes
mkdir -p $output-folder

for f in NC_017518_6k.fa ST154_6k.fa ST154Sim_6k.fa ST41Sim_6k.fa ST42Sim_6k.fa
do

    x=$(basename $f .fa)
    echo ${x}

    wgsim ${x}.fa -N 1000000 -1 150 -2 150  -e 0.005 -r 0 -R 0 -X 0 $output-folder/${x}.wgsim_er0.005.R1.fq  $output_folder/${x}.wgsim_er0.005.R2.fq
    gzip $output_folder/${x}.wgsim_er0.005.R1.fq
    gzip $output_folder/${x}.wgsim_er0.005.R2.fq

done
```

## build index for graph
```bash
mkdir refs

# copy graph to the refs work direvtory 
cp /home/zyang/pg_workshop/vg_deconstruct/4Sim_1K96.gfa /home/zyang/pg_workshop/graph_NGS/refs

#make tem_dir
mkdir /home/zyang/pg_workshop/graph_NGS/refs/temp_dir
```

```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      build_index_for_4SimGraph
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load vg/1.46.0

cd /home/zyang/pg_workshop/graph_NGS/refs
data=/home/zyang/pg_workshop/graph_NGS/refs/*.gfa
tem_dir=/home/zyang/pg_workshop/graph_NGS/refs/temp_dir

for f in $data
do

x=$(basename $f .gfa)

#convert the graph into 256 bp chunks, saving as vg format
vg mod -X 256 ${x}.gfa >${x}_256.vg

#build index of xg and gcsa index
vg index -b $tem_dir -t 48 -x ${x}_256.xg -g ${x}_256.gcsa -k 16 ${x}_256.vg

#small graph is ok without prunning, complex graph will need to prune first before generating index

### pruning: use -M if pruning fails
#vg prune -u -m node-mapping.tmp -t 48 -k 24 ${x}_256.vg > ${x}_256_chopped.vg

#vg index ${x}_256_chopped.vg -x ${x}_256_chopped.xg
### gcsa index
#vg index -b $tem_dir -t 48  -g ${x}_256_chopped.gcsa  ${x}_256_chopped.vg
done
```

## vg map NGS to graph 
```bash
mkdir /home/zyang/pg_workshop/graph_NGS/graph_based_mapping
```bash

```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      vgmap_5e_4Sim
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0


data=/home/zyang/pg_workshop/graph_NGS/simu_NGS_data/*R1.fq.gz
input_folder=/home/zyang/pg_workshop/graph_NGS/simu_NGS_data
output=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping


index=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.gcsa
basename=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256


for f in $data
do


x=$(basename $f R1.fq.gz)
echo ${x}

read1=${x}R1.fq.gz
read2=$(echo $read1|sed 's/R1.fq.gz/R2.fq.gz/')

echo $read2

#map paired reads using vg map
vg map -t 20  -d $basename -g $index  -f $input_folder/$read1 -f $input_folder/$read2 -N $x  > $output/${x}vgmap_4Sim.gam
#vg stats to check the mapping statistics
vg stats -a  $output/${x}vgmap_4Sim.gam  >$output/${x}vgmap_4Sim_stats

done
```


## genotying known variants 

```bash
mkdir /home/zyang/pg_workshop/graph_NGS/vgmap_12e_sim4_allR10S3_typing
```

generate snarls of graph
```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      vgmap_generate_snarls
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          1:00:00

module purge
module load vg/1.46.0


input_folder=/home/zyang/pg_workshop/graph_NGS/refs
output_folder=/home/zyang/pg_workshop/graph_NGS/refs

vg snarls $input_folder/4Sim_1K96_256.xg > $output_folder/4Sim_1K96_256.xg.snarls
```

genotyping 
```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      5e_vgmap_genotying
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0

data_gam=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping/*.wgsim_er0.005.vgmap_4Sim.gam
input=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping
output=/home/zyang/pg_workshop/graph_NGS/vgmap_12e_sim4_allR10S3_typing
graph_xg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg
snarls_file=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg.snarls



#compute snarls
#vg snarls $graph_xg >$snarls_file

for f in $data_gam
do

x=$(basename $f .wgsim_er0.005.vgmap_4Sim.gam)
echo ${x}


#Calculate the surpport reads ingoring mapping and base quality <5
#vg pack -t 48 -x $graph_xg -g $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -Q 5 -o $output/${x}vgmap_Sim4_256_aln.pack

#Calculate the surpport reads
vg pack -t 12 -x $graph_xg -g $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -o $output/${x}vgmap_sim4_256_aln.pack

#call variant using the same coordinates and including reference calls (for following compare)
vg call -t 12 -m 3,10 $graph_xg -k $output/${x}vgmap_sim4_256_aln.pack -r $snarls_file -a  >$output/${x}vgmap_sim4_256_aln.pack_allR10S3.vcf

done


```

## noverl variant calling using graph reference

```bash
mkdir /home/zyang/pg_workshop/graph_NGS/vgmap_5e_sim4_allR10S3_novelcalling 
```

```bash
#!/usr/bin/bash

#SBATCH --account       ga03793
#SBATCH --job-name      5e_vgmap_novelvariant_calling
#SBATCH --cpus-per-task 8
#SBATCH --mem           4G
#SBATCH --time          24:00:00

module purge
module load vg/1.46.0

data_gam=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping/*.wgsim_er0.005.vgmap_4Sim.gam
input=/home/zyang/pg_workshop/graph_NGS/graph_based_mapping
output=/home/zyang/pg_workshop/graph_NGS/vgmap_5e_sim4_allR10S3_novelcalling
graph_vg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.vg
graph_xg=/home/zyang/pg_workshop/graph_NGS/refs/4Sim_1K96_256.xg


#compute snarls
#vg snarls $graph_xg >$output/${graph_xg}.snarls

for f in $data_gam
do

x=$(basename $f .wgsim_er0.005.vgmap_4Sim.gam)
echo ${x}


#in order to also consider novel variants from the reads, use the augmented graph and gam (as created in the "Augmentation" example using vg augment -A)
#Augment augment the graph with all variation from the GAM, saving to aug.vg
### augment the graph with all variation from the GAM except
### that implied by soft clips, saving to aug.vg
### *aug-gam contains the same reads as aln.gam but mapped to aug.vg

vg augment -t 12 $graph_vg $input/${x}.wgsim_er0.005.vgmap_4Sim.gam -A $output/${x}nofilt_aug.gam >$output/${x}nofilt_aug.vg

#index the augmented graph
vg index -t 12 $output/${x}nofilt_aug.vg -x $output/${x}nofilt_aug.xg

## Compute the all read support from the augmented gam
vg pack -t 12 -x $output/${x}nofilt_aug.xg -g $output/${x}nofilt_aug.gam  -o $output/${x}nofilt_aug_allR.pack


#call variant
vg call -t 12 -m 3,10 $output/${x}nofilt_aug.xg -k $output/${x}nofilt_aug_allR.pack >$output/${x}nofilt_aug_allR.pack.vcf

#call variant snarl using the same coordinate
#vg call -t 48 -m 3,10 $output/${x}nofilt_aug.xg -k $output/${x}nofilt_aug_allR.pack -a >$output/${x}nofilt_aug_allR.pack_snarls.vcf

done
```



