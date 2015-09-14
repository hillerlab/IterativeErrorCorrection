# IterativeErrorCorrection

Iterative error correction of long 250 or 300 bp Illumina reads minimizes the total amount of erroneous reads, which improves contig assembly [1]. This pipeline runs multiple rounds of k-mer-based correction with an increasing k-mer size, followed by a final round of overlap-based correction. By combining the advantages of small and large k-mers, this pipeline is able to correct more base substitution errors, especially in repeats. The final overlap-based correction round can also correct small insertions and deletions. In [1], we show this higher read accuracy greatly improves contig assembly.

The script SGA-ICE (SGA-Iteratively Correcting Errors) implements iterative error correction by using modules from the String Graph Assembler (SGA) [2]. 

# Installation
First, you need to install this source code branch of SGA, where we added a few new parameters. 

 git clone https://github.com/ktrns/sga/commits/master

Then follow the [SGA installation instructions](https://github.com/jts/sga/tree/master/src#readme)

# Running SGA-ICE
All you need is a directory with the fastq files (ending *fastq or *fq). SGA-ICE creates a bash script with all commands for iterative error correction using default parameters that work well in general. To speedup the runtime, we recommended to set the number of threads to the number of cores available in your machine (-t num). 
 
 *SGA-ICE ***

Then, just executes the shell script. 
 *script.sh*
 
For each input file, you will get a final error corrected output file with the ending *finalec.fastq.

# Setting parameters
SGA-ICE allows to control all parameters if you do not want to use the default values. 

** Usage here **




# References
[1] Sameith K, Roscito J, Hiller M. Iterative error correction of long sequencing reads maximizes accuracy and improves contig assembly. Submitted

[2] Simpson JT and Durbin R (2012). Efficient de novo assembly of large genomes using compressed data structures. Genome Research, 22, 549–556.



