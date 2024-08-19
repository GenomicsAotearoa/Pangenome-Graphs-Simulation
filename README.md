# Pangenome Graphs - Simulation

This repository contains a section of the original [Pangenome Graphs Workshop](https://github.com/GenomicsAotearoa/Pangenome-Graphs-Workshop) that was removed to reduce the amount of content being covered.

It involves the use of the `wgsim` software to simulate sequencing a genome using short-read high throughput sequnecing technology. The short-read data can be used to investigate the difference in alignment percentages when aligning to a pangenome assembly, versus aligning to a single linear reference genome.


## Genomic data

The data for this exercise can be download by cloning this repository:

```
git clone https://github.com/mikblack/Pangenome-Graphs-Simulation
```

The data file `4Sim.fa` contains the genome sequence for four different genomes of *Neisseria meningitidis* (1 real, 3 with specific genomic feactures simulated), which need to be split into individual genomes. 

We can see the genome names via:

```
grep '^>' 4Sim.fa
```

The above command searches the file `5NM.fa` for lines that start ('^') with '>'.

We can split the file into the five separate genomes (using the '>' symbol to identify where each one starts) usign the following commands:

```
awk '/^>/ { gsub(">","",$1)
            FILE=$1 ".fa"
            print ">" $1 >> FILE
            next}
          { print >> FILE }' 4Sim.fa   
```

We now have four separate `.fa` (FASTA) files, one per genome:

```
NC_neisseria.fa
Sim1_3k.fa
Sim2_4k.fa
Sim3_5k.fa
```

## Use `wgsim` to simulate individual genomes and variation

Use `wgsim` to simulate 2x150bp NGS data, with an error rate of 0.005. 

```
output_folder=simulated_fastq
mkdir -p $output_folder

for f in NC_neisseria.fa Sim1_3k.fa Sim2_4k.fa Sim3_5k.fa
do
    x=$(basename $f .fa)
    echo ${x}

    wgsim -N 1000000 -1 150 -2 150  -e 0.005 -r 0 -R 0 -X 0 ${x}.fa $output_folder/${x}.wgsim_er0.005.R1.fq $output_folder/${x}.wgsim_er0.005.R2.fq
    gzip $output_folder/${x}.wgsim_er0.005.R1.fq
    gzip $output_folder/${x}.wgsim_er0.005.R2.fq

done
```

`wgsim` parameters:

 - `-N`: number of read pairs (1000000)
 - `-1`: length of the first read (150bp)
 - `-2`: length of the second read (150bp)
 - `-e`: per base sequencing error rate (0.005)
 - `-r`: de novo mutation rate (0)
 - `-R`: fractio of the genome comprising InDels (0)
 - `-X`: probability an indel is extended (0)

The simulated data from each individual can now be aligned to either a linear reference or a pangenome graph.

## Build index for graph

```
#convert the graph into 256 bp chunks, saving as vg format
vg mod -X 256 4Sim_1K96.gfa > 4Sim_1K96_256.vg
```
```
# Make temporary directory
mkdir TMP
#build index of xg and gcsa index
vg index -b TMP -t 48 -x 4Sim_1K96_256.xg -g 4Sim_1K96_256.gcsa -k 16 4Sim_1K96_256.vg
```

Small graph is ok without prunning, complex graph will need to prune first before generating index

Not run?

```
# pruning: use -M if pruning fails
# vg prune -u -m node-mapping.tmp -t 48 -k 24 ${x}_256.vg > ${x}_256_chopped.vg

# vg index ${x}_256_chopped.vg -x ${x}_256_chopped.xg
# gcsa index
# vg index -b $tem_dir -t 48  -g ${x}_256_chopped.gcsa  ${x}_256_chopped.vg
```

## vg map NGS to graph

```
data=simulated_fastq/*R1.fq.gz
input_folder=simulated_fastq
output=graph_based_mapping

mkdir -p $output

index=4Sim_1K96_256.gcsa
basename=4Sim_1K96_256

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

```
mkdir vgmap_12e_sim4_allR10S3_typing

# Generate snarls of graph
vg snarls 4Sim_1K96_256.xg > 4Sim_1K96_256.xg.snarls

# Genotyping
data_gam=graph_based_mapping/*.wgsim_er0.005.vgmap_4Sim.gam
input=graph_based_mapping
output=vgmap_12e_sim4_allR10S3_typing
graph_xg=4Sim_1K96_256.xg
snarls_file=4Sim_1K96_256.xg.snarls

#compute snarls
vg snarls $graph_xg >$snarls_file

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

## Novel variant calling using graph reference (some steps are the same as above)

```
mkdir vgmap_5e_sim4_allR10S3_novelcalling 

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
### *aug-gam contains the same reads as aln.gam but mapped to aug.vg

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

