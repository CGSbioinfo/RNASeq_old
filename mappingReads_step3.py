import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing
import subprocess

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raisec

def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths

###################
#  * Rscripts used:
# mapping_summary.R
# mapping_distribution.R
###################

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

logging.info(" ")
logging.info(" ")
logging.info("***************************************")
logging.info("*************MAPPING READS*************")
logging.info(" ")
logging.info("User command: " + str(sys.argv))


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... READING PARAMETERS FILE")
params_file = sys.argv[1]
args = {}
with open(params_file, 'r') as f:
    for line in f:
        entry = line.strip().split("=")
        if entry[0]:
            args[entry[0].strip(" ")] = entry[1].strip()
logging.info('-- Input parameters:')
path = args['Working directory'].replace("\\", "")
logging.info('Working Directory = ' + path)
gtfFile = args['GTF File']
logging.info('GTF File = ' + gtfFile)
refGenome = args['Reference Genome']
logging.info('Reference Genome = ' + refGenome)
nsamples = int(args['Number of samples'])
logging.info('Number of samples = ' + str(nsamples))
bedFile = args['BedFile']
logging.info('BedFile = ' + bedFile)
bedFile_10k = args['BedFile10K']
logging.info('BedFile_10k = ' + bedFile_10k)
refFlat = args['refFlat']
logging.info('refFlat = ' + refFlat)
rRNA_interval_list = args['rRNA_interval_list']
logging.info('rRNA_interval_list = ' + rRNA_interval_list)
strand = args['strand']
logging.info('strand = ' + strand)

os.chdir(path)

logging.info(" ")
logging.info(" ")
logging.info("#################################")
# Read in sampleNames
logging.info('... Sample names:')
sampleNames = []
sample_names_file = open('sample_names.txt','r')
for line in sample_names_file:
    sampleNames.append(line.strip())
    logging.info(line.strip())


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Mapping of reads")
version=subprocess.check_output('STAR --version', shell=True).strip()
logging.info("STAR version: "+ str(version))
make_sure_path_exists("./alignedReads/1pass_ann")
alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".sortedByCoord.out.bam", x)] # Check if aligned reads files exist

def alignment(i):
    trimmedReads = os.listdir("./trimmedReads")
    trimmedReads = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall(i, x)]
    r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2_.*val_2.fq", x)]
    if r2:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1_.*val_1.fq", x)]
        r2 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall("_R2_.*val_2.fq", x)]
        os.system("STAR --genomeDir " + refGenome + " --readFilesIn trimmedReads/" + r1[0] + " trimmedReads/" + r2[0] + " --runThreadN 4 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignedReads/1pass_ann/" + i)
        logging.info('Sample ' + i + ' done, pairedEnd mode')
    else:
        r1 = [trimmedReads[y] for y, x in enumerate(trimmedReads) if re.findall( "_R1_.*trimmed.fq", x)]
        os.system("STAR --genomeDir " + refGenome + " --readFilesIn trimmedReads/" + r1[0] + " --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix alignedReads/1pass_ann/" + i)
        logging.info('Sample ' + i + ' done, singleEnd mode')
if not indicesAlignmentFiles:
    Parallel(n_jobs=2)(delayed(alignment)(i) for i in sampleNames)
logging.info("... Finished Mapping")

# Index Bam files
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Indexing BAM files")
version=subprocess.check_output(' samtools --version | grep samtools', shell=True).strip()
logging.info("Samtools version: "+ str(version))
alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".sortedByCoord.out.bam.bai", x)] # Check if aligned reads files exist
def indexing(i):
    os.system("samtools index alignedReads/1pass_ann/" + i + "Aligned.sortedByCoord.out.bam")
    logging.info('Sample ' + i + ' done')
if not indicesAlignmentFiles:
    Parallel(n_jobs=7)(delayed(indexing)(i) for i in sampleNames)
logging.info("... Finished Indexing")


# Mapping QC
logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Mapping QC")
os.system("mapping_summary.R .")
logging.info('Mapping summary done')
logging.info(" ")

# GeneBody Coverage
make_sure_path_exists("./alignedReads/1pass_ann/QC")
make_sure_path_exists("./Report/figure/")
alignmentQCOut = os.listdir(path+"/alignedReads/1pass_ann/QC") # Check if aligned reads files exist
indicesAlignmentQCFiles = [i for i, x in enumerate(alignmentQCOut) if re.findall(".geneBodyCoverage.", x)] # Check if aligned reads files exist
if not indicesAlignmentQCFiles:
    os.system("ls ./alignedReads/1pass_ann/*Aligned.sortedByCoord.out.bam > tempbamfiles.txt")
    os.system("python ~/bin/geneBody_coverage.py -r " + bedFile_10k + " -i tempbamfiles.txt -o alignedReads/1pass_ann/QC/10KGenes")
    os.system("rm tempbamfiles.txt")
# os.system('cp alignedReads/1pass_ann/QC/10KGenes.geneBodyCoverage.curves.pdf Report/figure/10KGenes_geneBodyCoverage_curves.pdf')
logging.info('GeneBody Coverage done')
logging.info(" ")

# Junctions and junction saturation
def junctions(i):
    os.system("python ~/bin/junction_annotation.py -i alignedReads/1pass_ann/" + i +"Aligned.sortedByCoord.out.bam  -o alignedReads/1pass_ann/QC/" + i + " -r " + bedFile + " & ")
    os.system("python ~/bin/junction_saturation.py -i alignedReads/1pass_ann/" + i +"Aligned.sortedByCoord.out.bam  -o alignedReads/1pass_ann/QC/" + i + " -r " + bedFile)
    logging.info("Calculated junction annotation and junction saturation sample " + i)
qcOut = os.listdir(path+"/alignedReads/1pass_ann/QC")
qcFiles = [i for i, x in enumerate(qcOut) if re.findall(".junction.", x)] # Check if aligned reads files exist
if not qcFiles:
     Parallel(n_jobs=8)(delayed(junctions)(i) for i in sampleNames)
os.system('grep "y=c(" alignedReads/1pass_ann/QC/*junctionSaturation*  | sed \'s/:y=c(/,/g\' | sed \'s/.junctionSaturation_plot.r//g\' | sed \'s/)//g\' | sed \"s/.*\///g\"  > alignedReads/1pass_ann/QC/junctionSat_all.csv')
os.system('~/bin/junctionPlotAll.R .')
logging.info(" ")

# Collect Metrics
qcOut = os.listdir(path+"/alignedReads/1pass_ann/QC") # Check if aligned reads files exist
qcFiles = [i for i, x in enumerate(qcOut) if re.findall(".metrics.txt", x)] # Check if aligned reads files exist
logging.info("Calculating RNASeq Metrics")
#version=subprocess.check_output("java -jar ~/tools/picard-tools-1.127/picard.jar CollectRnaSeqMetrics --version",shell=TRUE).strip()
#version=version.split("(")[0]
logging.info("Samtools version: "+ str(version))
if not qcFiles:
    os.system("find alignedReads/1pass_ann/*.sortedByCoord.out.bam | sed \"s/.*\///g\" | sed \"s/.sortedByCoord.out.bam//g\" | " # get file names and format them
          "parallel -j 3 --no-notice "
          "\"java -jar ~/tools/picard-tools-1.127/picard.jar CollectRnaSeqMetrics "
          "REF_FLAT=" + refFlat + " "
          "RIBOSOMAL_INTERVALS=" + rRNA_interval_list + " "
          "STRAND_SPECIFICITY=" + strand + " "
          "INPUT=alignedReads/1pass_ann/{}.sortedByCoord.out.bam "
          "OUTPUT=alignedReads/1pass_ann/QC/{}_metrics.txt \"")
    logging.info("Calculated RNASeq Metrics")
def pct(i):
    os.system('mapping_distribution.R ./alignedReads/1pass_ann/QC/ ' + i)
Parallel(n_jobs=7)(delayed(pct)(i) for i in sampleNames)
logging.info("... Finished mapping QC")


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Sorting BAM files by name ")
version=subprocess.check_output(' samtools --version | grep samtools', shell=True).strip()
logging.info("Samtools version: "+ str(version))
# Sort bam files by name
alignmentOut = os.listdir(path+"/alignedReads/1pass_ann") # Check if aligned reads files exist
indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".sortedByCoord.sortedByName.out.bam", x)] # Check if aligned reads files exist
def sort(i):
    os.system('samtools sort -n alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.out.bam alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.sortedByName.out')
    logging.info("Sample " + i + ' done')
if not indicesAlignmentFiles:
    Parallel(n_jobs=8)(delayed(sort)(i) for i in sampleNames)
logging.info("... Finished sorting BAM files by name")

logging.info(" ")
logging.info(" ")
logging.info(" ")
logging.info("##################################################################")
logging.info("------------------------------------------------------------------")
logging.info("##################################################################")
logging.info(" ")
logging.info(" ")
logging.info(" ")

