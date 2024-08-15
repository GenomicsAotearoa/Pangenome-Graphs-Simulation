# Pangenome Graphs - Simulation

This repository contains a section of the original [Pangenome Graphs Workshop](https://github.com/GenomicsAotearoa/Pangenome-Graphs-Workshop) that was removed to reduce the amount of content being covered.

It involves the use of the `wgsim` software to simulate sequencing a genome using short-read high throughput sequnecing technology. The short-read data can be used to investigate the difference in alignment percentages when aligning to a pangenome assembly, versus aligning to a single linear reference genome.


## Genomic data

The data for this exercise can be download by cloning this repository:

```
git clone https://github.com/ZoeYang2020/dataset_for_pg_workshop
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
output_folder=simulated_fastq
mkdir -p $output_folder

for f in NC_017518.fa ST154Sim.fa ST41Sim.fa ST42Sim.fa
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

The simulated data from each individual can now be aligned to either a linear reference or a pangenome graph using the methods detailed in 
