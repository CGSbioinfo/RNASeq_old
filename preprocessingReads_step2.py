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
# preprocessing_numbers.R
# trimgalore_summary.R
# mapping_summary.R
###################

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')


logging.info("*************PREPROCESSING-READS*************")
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

# Check if trimmed files exist, might be removed
trimmedFiles = []
for root, dir, files in os.walk('./trimmedReads/'):
   trimmedFiles.extend(files)
indicesTrimmedFiles = [trimmedFiles[i] for i, x in enumerate(trimmedFiles) if re.findall(".fastq.gz", x)]

# Check if input files are uncompressed or .gz, at this point all files should be compressed and NOT SPLITTED BY LANE
readFiles = []
for root, dir, files in os.walk(path):
    readFiles.extend(files)
indicesgzFiles = [i for i, x in enumerate(readFiles) if re.findall(".fastq.gz", x)]
logging.info(' ')
if indicesgzFiles:
   gz =".gz"
   logging.info("Input files are .gz ")
else:
   gz =""
   logging.info("Input files are uncompressed ")

logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... PRE-PROCESSING OF READS")
version=subprocess.check_output('fastqc --version', shell=True).strip()
logging.info("Fastqc version: "+ str(version))
version=subprocess.check_output('trim_galore --version | grep "version"', shell=True).strip()
logging.info("TrimGalore version: "+ str(version))
version=subprocess.check_output('cutadapt --version ', shell=True).strip()
logging.info("cutadapt version: "+ str(version))
make_sure_path_exists("./trimmedReads")
def preprocessing(i):
    allFiles = os.listdir("./rawReads/" + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2_", x)]
    os.system("fastqc ./rawReads/" + i + "/" + i + "_R1_*.fastq" + gz + " --outdir=rawReads/" + i)
    os.system("unzip ./rawReads/" + i + "/" + i + "_R1_*_fastqc.zip -d rawReads/" + i + "/")
    if pairedReads_temp:
        os.system("fastqc ./rawReads/" + i + "/" + i + "_R2_*.fastq" + gz + " --outdir=rawReads/" + i)
        os.system("unzip ./rawReads/" + i + "/" + i + "_R2_*_fastqc.zip -d rawReads/" + i + "/")
        os.system("trim_galore --paired --retain_unpaired --fastqc --output_dir trimmedReads/ " + 'rawReads/' + i + "/" + i + "_R1_*.fastq" + gz + " rawReads/" + i + "/" + i + "_R2_*.fastq" + gz)
        os.system("unzip trimmedReads/" + i + "_R1_*val_1_fastqc.zip -d trimmedReads/")
        os.system("unzip trimmedReads/" + i + "_R2_*val_2_fastqc.zip -d trimmedReads/")
    else:
        os.system("trim_galore --fastqc --output_dir trimmedReads/ rawReads/" + i + "/" + i + "_R1_*.fastq" + gz)
        os.system("unzip trimmedReads/" + i + "_R1_*_trimmed_fastqc.zip -d trimmedReads/")
    logging.info("Preprocessing done for sample " + i)

# Running jobs
if not indicesTrimmedFiles:
	Parallel(n_jobs=7)(delayed(preprocessing)(i) for i in sampleNames)

logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Creating files with output information and summary")
# Create files with output information (used in the Report)
os.system("preprocessing_numbers.R .")
os.system("trimgalore_summary.R .")

# Create files from fastqc_data.txt for creating plots (used in the Report)
def tables(i):
    os.system('fastqc_plot_data.py ./rawReads/' + i + '/' + i + '_R1_*fastqc/fastqc_data.txt all ./rawReads/' + i + '/' + i + '_R1_*fastqc/')
    allFiles = os.listdir("./rawReads/" + i )
    pairedReads_temp = [allFiles[y] for y, x in enumerate(allFiles) if re.findall("_R2_", x)]
    if pairedReads_temp:
        os.system('fastqc_plot_data.py ./rawReads/' + i + '/' + i + '_R2_*fastqc/fastqc_data.txt all ./rawReads/' + i + '/' + i + '_R2_*fastqc/')
        os.system('fastqc_plot_data.py ./trimmedReads/' + i + '_R1_*val_1_fastqc/fastqc_data.txt all ./trimmedReads/' + i + '_R1_*val_1_fastqc/')
        os.system('fastqc_plot_data.py ./trimmedReads/' + i + '_R2_*val_2_fastqc/fastqc_data.txt all ./trimmedReads/' + i + '_R2_*val_2_fastqc/')
    else:
        os.system('fastqc_plot_data.py ./trimmedReads/' + i + '_R1_*trimmed_fastqc/fastqc_data.txt all ./trimmedReads/' + i + '_R1_*trimmed_fastqc/')
Parallel(n_jobs=8)(delayed(tables)(i) for i in sampleNames)

os.system('ls rawReads/*/*fastqc  |  grep -v trimmed  | grep ":"  | sed \'s/://g\' > sample_names2.txt')
os.system('fastqc_summary.py ./sample_names2.txt ./summary_fastqc.txt')

logging.info(" ")
logging.info(" ")
logging.info(" ")
logging.info("##################################################################")
logging.info("------------------------------------------------------------------")
logging.info("##################################################################")
logging.info(" ")
logging.info(" ")
logging.info(" ")
#logging.info("DONE PRE-PROCESSING OF READS")

