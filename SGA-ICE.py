#!/usr/bin/env python

""" SGA-ICE
SGA Iteratively Correcting Errors
"""
__author__ = ["Juliana Roscito, Katrin Sameith, Michael Hiller"] 

import os
import stat
import sys
import argparse

def __parse_arguments():
    """Read arguments from the command line.
    Returns an argument object containing all input parameters.
    """

    parser = argparse.ArgumentParser(description='SGA-ICE produces a shell script \'runMe.sh\' that contains all commands to run iterative error correction of the given read data with the given parameters. \nRead data must be in fastq format and files need to have the ending .fastq or .fq.')
    parser.add_argument('inputDir', type=str, help='Path to directory with the *.fastq or *.fq files. The produced \'runMe.sh\' will be located here.')
    parser.add_argument('-k', '--kmers', type=str, default=None, help='List of k-mers for k-mer correction; values should be comma-separated. If -k is not provided, SGA-ICE does 3 rounds of k-mer correction with k-mer sizes determined based on the length of the read from the first file in inputDir. We advise the user to choose k-mer values manually if the sequences in the *.fastq files have different read lengths.')
    parser.add_argument('-t', '--threads', type=int, default='1', help='Number of threads used. Default is 1. Set to higher values if you have more than one core and want to reduce the runtime.')
    parser.add_argument('--noOvlCorr', action="store_true", help='If set, do not run a final overlap-based correction round.')
    parser.add_argument('--noCleanup', action="store_true", help='If set, keep all intermediate files in the temporary directory.')

    # print help if no argument given
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(2)

    args = parser.parse_args()

    print '### Input directory is:', args.inputDir
    if not os.path.isabs(args.inputDir):
        args.inputDir = os.path.join(os.getcwd(), args.inputDir)

    print '### Number of threads is:', args.threads

    if args.kmers:
        args.kmers = [int(x) for x in args.kmers.split(',')]
        print '### User provided k-mers are:', args.kmers
    else:
        print '### User did not provide k-mers'

    if args.noOvlCorr:
        print '### User chose not to do overlap correction'
    else:
        print '### User chose to do overlap correction'

    if args.noCleanup:
         print '### Temporary files will not be deleted'
    else:
         print '### All temporary files will be deleted'

    return args


def get_file_list(inputDir):
    """Get list of *.fastq and *.fq files in the input directory.
    Returns triple of lists of equal length:
        * list of *.fastq and *.fq files in the input directory
        * list of file prefixes
        * list of file endings
    """

    files = []
    filePrefix = []
    fileEnding = []
    for f in os.listdir(inputDir):
        if f.endswith('.fastq'):
            files.append(f)
            fprefix = f.split('.fastq')[0]
            filePrefix.append(fprefix)
            fileEnding.append("fastq")
        if f.endswith('.fq'):
            files.append(f)
            fprefix = f.split('.fq')[0]
            filePrefix.append(fprefix)
            fileEnding.append("fq")

    if not files:
        print "ERROR: There is no *.fastq or *.fq file in this directory: %s" % inputDir
        sys.exit(1) 

    print "### These are the fastq files found in " + inputDir + ":\n" + "\n".join(files)
    return files, filePrefix, fileEnding


def get_read_length(inputDir, exampleFile):
    """Get read length from one of the input files.
    Returns read length. 
    """

    with open(inputDir + "/" + exampleFile) as f:
        seq = f.readlines()[1].rstrip("\n")
        seqLength = len(seq)

    print ("### Your reads are %(seqLength)d bp long, calculated from file: %(exampleFile)s"
          % dict(seqLength=seqLength, exampleFile=os.path.join(inputDir,exampleFile)))
    return seqLength


def get_kmers(seqLength, kmers=None):
    """Get k-mer values for k-mer-correction. 
    Returns a list of k-mers. If k-mers are not provided, default is set to three ks. 
    """

    if not kmers:
        kmer1 = 40
        kmer3 = (seqLength*2)/3
        kmer2 = ((kmer3-kmer1)/2)+kmer1
        kmers = [kmer1, kmer2, kmer3]
        print '### Using default k-mer values for k-mer error correction:', kmers

    kmers.sort()
    # sanity checks: k-mers should be >0 and <= seqLength, no k-mer should occur twice
    if kmers[-1] > seqLength:
        print "ERROR: This k-mer value is greater than sequence length: %d" % kmers[-1]
        sys.exit(1)
    if kmers[0] <= 0:
        print "ERROR: This k-mer value is less or equal than zero: %d" % kmers[0]
        sys.exit(1)
    for i in range(1, len(kmers)):
        if kmers[i] == kmers[i-1]:
            print "ERROR: Two k-mers in the list are equal: %d" % kmers[i]
            sys.exit(1)
    return kmers


def sga_ice_write(args, filePrefix, fileEnding):
    """Write output 'runMe.sh' script. This is the script that iteratively runs 
    SGA modules. 
    Returns True when writing finished successfully.
    """ 


    # add 0 to the list of k-mers. The input into the first round will be called *.k0.fastq
    kmers = [0] + args.kmers
    overlapCorrection = not args.noOvlCorr
    cleanup = not args.noCleanup

    values = {"inputDir": args.inputDir,
              "threads": args.threads,
              "kmers": kmers,
              "kmers[-1]": kmers[-1]
             }

    # create runMe.sh script
    sgaICE = open(args.inputDir + "/runMe.sh", "w")
    sgaICE.write(("#!/bin/bash\n"
                  "set -e\n" 
                  "set -o pipefail\n\n\n"
                  "### Create 'ec' directory for the final error-corrected files ###\n"
                  "mkdir -p %(inputDir)s/ec\n\n"
                  "### Create a unique temporary directory to keep the temporary files (preprocessed files, files after each correction round) ###\n"
                  "tmpDir=`mktemp -d`; cd $tmpDir\n"
                  "echo \"Temporary directory is: $tmpDir\"\n"
                  "\n### Link input files to temp dir ###\n"
                  "ln -s %(inputDir)s/*.f*q .\n\n")
                 % values)


    # add commands for SGA-preprocessing each input file
    sgaICE.write("### Run preprocessing for each input file ###\n"
                 "echo \'### Start sga preprocessing ###\'\n")
    for (file, ending) in zip(filePrefix, fileEnding): 
        sgaICE.write(("sga preprocess "
                      "--no-primer-check "
                      "--pe-mode 0 "
                      "--permute-ambiguous "
                      "--min-length 0"
                      "--out %(file)s.pp.%(ending)s "
                      "%(file)s.%(ending)s \n")
                     % dict(values, file=file, ending=ending))

    # SGA-index the preprocessed input data
    sgaICE.write("echo \'    Build new index from the preprocessed input files\'\n")
    # link it if we have a single input file, otherwise cat them all
    if len(filePrefix) == 1:
        sgaICE.write("ln -s *.pp.f*q all.pp.fastq\n")
    else:
        sgaICE.write("cat *.pp.f*q > all.pp.fastq\n")

    sgaICE.write("sga index -a ropebwt -t %(threads)d all.pp.fastq\n"
                  % values)

    sgaICE.write("rename .pp. .ec.k0. *.pp.*\n"
                 "echo \'### Finished sga preprocessing ###\'\n\n")


    # output all commands for all k-mer correction rounds
    sgaICE.write(("### Do k-mer based correction\n"
                 "echo \'### Starting %s rounds of k-mer based error correction with the following ks: %s  ###\'\n")
                 % (len(kmers)-1, ", ".join([str(k) for k in kmers[1:]])))

    for kIdx in range(1, len(kmers)):
        sgaICE.write("echo \'  Running correction round with k=%d \'\n" % kmers[kIdx])
        for (file, ending) in zip(filePrefix, fileEnding): 
            prevIdx = kIdx-1
            sgaICE.write(("echo \'    Correcting file %(file)s.ec.k%(prevK)d.%(ending)s\'\n"
                          "sga correct "
                          "--count-offset 2 "
                          "--threads %(threads)d "
                          "--prefix all.ec.k%(prevK)d "
                          "--outfile %(file)s.ec.k%(k)d.%(ending)s "
                          "-k %(k)d "
                          "--learn "
                          "%(file)s.ec.k%(prevK)d.%(ending)s\n")
                         % dict(values,
                                file=file,
                                ending=ending,
                                prevK=kmers[prevIdx],
                                k=kmers[kIdx]))


        # index the newly corrected reads, unless it is the last round of k correction
        if kIdx == kmers.index(kmers[-1]):
            sgaICE.write("echo \'### Finished all k-mer based correction rounds ###\'\n\n")
        else:
            sgaICE.write("echo \'    Build new index\'\n")
            # link it if we have a single input file, otherwise cat them all
            if len(filePrefix) == 1:
                sgaICE.write(("ln -s *.ec.k%(k)d.f*q all.ec.k%(k)d.fastq\n")
                              % dict(k=kmers[kIdx]))
            else:
                sgaICE.write(("cat *.ec.k%(k)d.f*q > all.ec.k%(k)d.fastq\n")
                              % dict(k=kmers[kIdx]))
            sgaICE.write(("sga index -a ropebwt -t %(threads)d all.ec.k%(k)d.fastq\n")
                          % dict(values, k=kmers[kIdx]))


    # output commands for the final overlap-based correction round
    if overlapCorrection:
        sgaICE.write("### Do overlap-based correction\n"
                      "echo \'### Run overlap-based correction ###\'\n"
                      "echo \'    Build new index from the output of k-mer correction\'\n")
        # link it if we have a single input file, otherwise cat them all
        if len(filePrefix) == 1:
            sgaICE.write(("ln -s *.ec.k%(kmers[-1])d.f*q all.ec.k%(kmers[-1])d.fastq\n")
                         % values)
        else:
            sgaICE.write(("cat *.ec.k%(kmers[-1])d.f*q > all.ec.k%(kmers[-1])d.fastq\n")
                         % values)
        sgaICE.write((	"sga index -a ropebwt -t %(threads)d all.ec.k%(kmers[-1])d.fastq\n")
                      % values)

        for (file,ending) in zip(filePrefix, fileEnding): 
            sgaICE.write(("echo \'    Overlap-correction of file %(file)s.ec.k%(kmers[-1])d.%(ending)s\'\n"
                          "sga correct "
                          "--prefix all.ec.k%(kmers[-1])d "
                          "--algorithm overlap "
                          "--error-rate 0.01 "
                          "-m 40 --threads %(threads)d "
                          "--outfile %(file)s.final.ecOv.%(ending)s "
                          "%(file)s.ec.k%(kmers[-1])d.%(ending)s\n")
                         % dict(values, file=file, ending=ending))
        sgaICE.write("echo \'### Finished overlap-based correction ###\'\n\n")

        sgaICE.write("### After overlap correction, .fastq files are converted to .fasta because the quality values in the .fastq file produced by 'sga correct' not always match the length of the reads if insertions/deletions were corrected.\n"
                     "### Newer versions of sga may solve this issue\n"
                     "echo \'### After overlap correction, convert fastq to fasta files ###\'\n")

        for (file,ending) in zip(filePrefix, fileEnding): 
            # convert fastq to fasta:
            sgaICE.write(("awk 'BEGIN {OFS = \"\\n\"} "
                          "{header = $0 ; getline seq ; getline qheader ; getline qseq ; print \">\"header, seq}' "
                          "< %(file)s.final.ecOv.%(ending)s > %(file)s.final.ecOv.%(ending)s.fasta\n")
                         % dict(values, file=file, ending=ending))

        # move the final files afte overlap-correction
        sgaICE.write(("mv *final.ecOv.f*q* %(inputDir)s/ec\n"
                      "\necho; echo; echo \'### ALL DONE.  Final error-corrected files after k-mer-based and overlap-based correction are %(inputDir)s/ec/*.final.ecOv.f*q* ###\'\n")
                     % values)

    else:
        # move the final files after k-mer correction
        sgaICE.write("echo '### Copying final files to \'ec\' directory ###\'\n")
        for (file,ending) in zip(filePrefix, fileEnding): 
            sgaICE.write(("cp %(file)s.ec.k%(kmers[-1])d.%(ending)s  %(inputDir)s/ec/%(file)s.final.ecKmer.%(ending)s\n")
                         % dict(values, file=file, ending=ending))
        sgaICE.write(("\necho; echo; echo '### ALL DONE.  Final error-corrected files after k-mer-based correction are %(inputDir)s/ec/*.final.ecKmer.f*q ###\'\n")
                     % values)

    if cleanup:
        sgaICE.write("\n### Remove temp dir: ###\n"
                     "rm -rf $tmpdir\n")
    else:
        sgaICE.write("echo \"### SGA-ICE is finished. All temporary files are in this directory: $tmpDir\"\n")

    sgaICE.close()




    return True



if __name__ == "__main__":
    args = __parse_arguments()
    files, filePrefix, fileEnding = get_file_list(args.inputDir)

    seqLength = get_read_length(args.inputDir, files[0])
    args.kmers = get_kmers(seqLength, args.kmers)

    sga_ice_write(args, filePrefix, fileEnding)

    # make runMe.sh executable
    st = os.stat(args.inputDir + "/runMe.sh")
    os.chmod(args.inputDir + "/runMe.sh", st.st_mode | stat.S_IEXEC)

    print "\n\nTo run %d correction round with k=%s" % (len(args.kmers), ",".join([str(k) for k in args.kmers])),
    if not args.noOvlCorr:
        print "and a final round of overlap-based correction",
    print "using %d threads, execute \n%s/runMe.sh" % (args.threads, args.inputDir)
