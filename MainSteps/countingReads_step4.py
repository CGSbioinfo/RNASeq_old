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
###################

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')


logging.info("*************Counting Reads*************")
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
#logging.info("\n -- Input parameters: \n Working Directory = " + path + "\n GTF file = " + gtfFile + "\n Reference Genome = " + refGenome + "\n Number of samples = " + str(nsamples) + "\n BedFile = " + bedFile + "\n BedFile_10k = " + bedFile_10k + "\n refFlat = " + refFlat + "\n rRNA_interval_list = " + rRNA_interval_list + "\n strand = " + strand )

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
logging.info("..Counting reads")
version=subprocess.check_output('htseq-count -h | grep version | sed \'s/.*,//g\' | sed \'s/^ //g\' ', shell=True).strip()
logging.info("HTSeq version: "+ str(version))
make_sure_path_exists("./countedReads/1pass_ann/")
countOut = os.listdir(path+"/countedReads/1pass_ann/")
indicesCountFiles = [i for i, x in enumerate(countOut) if re.findall(".count", x)]

if strand == 'SECOND_READ_TRANSCRIPTION_STRAND':
    s_htseq = 'reverse'
else:
    s_htseq = 'yes'

def count(i):
    os.system('samtools view alignedReads/1pass_ann/' + i + 'Aligned.sortedByCoord.sortedByName.out.bam | htseq-count -a 10 -m union -s ' + s_htseq + ' - ' + gtfFile + ' > countedReads/1pass_ann/' + i + '.count')
    logging.info('Sample ' + i + ' done, ' + s_htseq + ' mode')

if not indicesCountFiles:
    Parallel(n_jobs=8)(delayed(count)(i) for i in sampleNames)

os.system("~/bin/countsLog_rnaseq.R .")

logging.info("... Finished Counting reads")

logging.info(" ")
logging.info(" ")
logging.info(" ")
logging.info("##################################################################")
logging.info("------------------------------------------------------------------")
logging.info("##################################################################")
logging.info(" ")
logging.info(" ")
logging.info(" ")
