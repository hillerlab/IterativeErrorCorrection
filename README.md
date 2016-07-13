# IterativeErrorCorrection

Iterative error correction of long 250 or 300 bp Illumina reads minimizes the total amount of erroneous reads, which improves contig assembly [1]. This pipeline runs multiple rounds of k-mer-based correction with an increasing k-mer size, followed by a final round of overlap-based correction. By combining the advantages of small and large k-mers, this pipeline is able to correct more base substitution errors, especially in repeats. The final overlap-based correction round can also correct small insertions and deletions. In [1], we show this higher read accuracy greatly improves contig assembly.

The script SGA-ICE (SGA-Iteratively Correcting Errors) implements iterative error correction by using modules from the String Graph Assembler (SGA) [2]. 

# Installation
First, you need to install SGA version v0.10.14 or later. 

 `git clone https://github.com/jts/sga.git`

Then follow the [SGA installation instructions](https://github.com/jts/sga/tree/master/src#readme)

# Running SGA-ICE
All you need is a directory with the fastq files (ending *fastq or *fq). SGA-ICE creates a 'runMe.sh' script with all commands for iterative error correction using default parameters that work well in general. To speedup the runtime, we recommended to set the number of threads to the number of cores available in your machine (-t num). 

**Example:** 
  `SGA-ICE.py /path/to/fastq/data/ -t 8`

If you are happy with the default parameters, just execute 'runMe.sh'. 
 
The error corrected files will be located in a /path/to/fastq/data/ec directory.

# Setting parameters
SGA-ICE allows to control all parameters if you do not want to use the default values. 

```
usage: SGA-ICE.py [-h] [-k KMERS] [-t THREADS] [--noOvlCorr] [--noCleanup]
                  [--scriptName SCRIPTNAME] [--errorRate ERRORRATE]
                  [--minOverlap MINOVERLAP]
                  inputDir

SGA-ICE produces a shell script that contains all commands to run iterative
error correction of the given read data with the given parameters. Read data
must be in fastq format and files need to have the ending .fastq or .fq.

positional arguments:
  inputDir              Path to directory with the *.fastq or *.fq files. The
                        produced shell script will be located here.

optional arguments:
  -h, --help            show this help message and exit
  -k KMERS, --kmers KMERS
                        List of k-mers for k-mer correction; values should be
                        comma-separated. If -k is not provided, SGA-ICE does 3
                        rounds of k-mer correction with k-mer sizes determined
                        based on the length of the read from the first file in
                        inputDir. We advise the user to choose k-mer values
                        manually if the sequences in the *.fastq files have
                        different read lengths.
  -t THREADS, --threads THREADS
                        Number of threads used. Default is 1. Set to higher
                        values if you have more than one core and want to
                        reduce the runtime.
  --noOvlCorr           If set, do not run a final overlap-based correction
                        round.
  --noCleanup           If set, keep all intermediate files in the temporary
                        directory.
  --scriptName SCRIPTNAME
                        Name of the shell script containing the error
                        correction commands. By default, script is called
                        runMe.sh
  --errorRate ERRORRATE
                        sga correct -e parameter for overlap correction.
                        Maximum error rate allowed between two sequences to
                        consider them overlapped. Default is 0.01
  --minOverlap MINOVERLAP
                        sga correct -m parameter for overlap correction.
                        Minimum overlap required between two reads. Default is
                        40
```

**Example:**
  `SGA-ICE.py /path/to/fastq/data/ -k 40,60,100,125,150,200 --noCleanup --noOvlCorr --scriptName correctMyData.sh`



# References
[1] Sameith K, Roscito J, Hiller M (2016). [Iterative error correction of long sequencing reads maximizes accuracy and improves contig assembly](http://bib.oxfordjournals.org/content/early/2016/02/10/bib.bbw003.long). Briefings in Bioinformatics, doi: 10.1093/bib/bbw003

[2] Simpson JT and Durbin R (2012). Efficient de novo assembly of large genomes using compressed data structures. Genome Research, 22, 549-556.


# Comments, Requests, Bug reports
Please email `hiller@mpi-cbg.de`
