# goal is to create distant paired ends using my paired end finder
# then use output to generate far junction library
# then use output to align to bowtie

# input orig files (can't have other files!)
# name the data set.  the output file will be created with this name and some underscores.
# name output directory
# minimum number of base pairs apart that user wants to find paired ends
# need pickle file for

def checkProcesses(popenDict):
	"""
	Function :
	Args     : popenDict - A dict. whose keys are subprocess.Popen instances. The value of a key is also a dict. that has the following
								keys: 'stdout', 'stderr', and 'cmd'.
	"""
	for popen in popenDict:
		popen.communicate() #hangs until job finishes
		retcode = popen.returncode
		if retcode:
			cmd = popenDict[popen]['cmd']
			stdout = popenDict[popen]['stdout']
			stderr = popenDict[popen]['stderr']
			raise Exception("Command '{cmd}' failed with return code {retcode}. Log files are {stdout} and {stderr}.".format(cmd=cmd,stdout=stdout,stderr=stderr))

import os
import sys
import subprocess
from argparse import ArgumentParser
description = ""

parser = ArgumentParser(description=description)
parser.add_argument("CIRCPIPE_DIR",help="Dir containing circ pipe output (incl linda's directories orig, circReads, logs, sample stats)")
parser.add_argument("OUTPUT_DIR",help="Far Junc Dir")
parser.add_argument("USERBPDIST",help="using 100000 (100Kb) so far for testing purposes")
parser.add_argument("REFGENOME",help="HG19 vs HG38;could upgrade to HG38.")
parser.add_argument("NUMBASESAROUNDJUNC",help="default for linda's is 8 for read lengths < 70 and 13 for read lengths > 70")
parser.add_argument("NumIndels",help="current default = 5, arbitrary")

args = parser.parse_args()

CIRCPIPE_DIR = args.CIRCPIPE_DIR
OUTPUT_DIR = args.OUTPUT_DIR
USERBPDIST = 100000 #args.USERBPDIST
REFGENOME = "HG19" #args.REFGENOME
NUMBASESAROUNDJUNC = args.NUMBASESAROUNDJUNC
NumIndels = 5 #args.NumIndels

#end arg parsing


## REPLACE THESE THREE FIELDS AFTER INSTALLATION
MACHETE="/src/machete" #nathankw - formerly called INSTALLDIR
CIRCREF="/share/PI/horence/circularRNApipeline_Cluster/index" #nathankw - update this to path to reference libraries output by KNIFE (directory that contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and hg19_junctions_scrambled bowtie indices). Probably will need to set this as a runtime parameter.
REG_INDEL_INDICES="/home/data/IndelIndices"


if REFGENOME == "HG38":
	EXONS="/scratch/PI/horence/grch38_junctions"
elif REFGENOME == "HG19":
	EXONS="/home/data/HG19exons" #nathankw - formerly called PICKLEDIR
else:
	raise Exception("Incorrect value for REFGENOME. Must be one of HG19 or HG38.")


ORIG_DIR=os.path.join(CIRCPIPE_DIR,"orig")
GLM_DIR=os.path.join(CIRCPIPE_DIR,"circReads/glmReports")
DistantPEDir=os.mkdir(os.path.join(OUTPUT_DIR,"DistantPEFiles"))
FASTADIR=os.mkdir(os.path.join(OUTPUT_DIR,"fasta"))
BOWTIE_DIR=os.mkdir(os.path.join(OUTPUT_DIR,"BowtieIndex"))
UNALIGNEDDIR=os.path.join(ORIG_DIR),"unaligned")
FARJUNCDIR=os.mkdir(os.path.join(OUTPUT_DIR,"FarJunctionAlignments"))
SECONDFARJUNCDIR=os.mkdir(os.path.join(OUTPUT_DIR,"FarJuncSecondary"))
BadFJDir=os.mkdir(os.path.join(OUTPUT_DIR,"BadFJ"))
StemFile=os.path.join(OUTPUT_DIR,"StemList.txt")


os.mkdir(os.path.join(OUTPUT_DIR,"reports"))
LOG_DIR = os.mkdir(os.path.join(OUTPUT_DIR,"err_and_out"))
os.mkdir(os.path.join(OUTPUT_DIR,"BowtieIndels"))
os.mkdir(os.path.join(OUTPUT_DIR,"FarJuncIndels"))
os.makedirs(os.path.join(SECONDFARJUNCDIR,"AlignedIndels/RemoveNonOverlap"))
os.mkdir(os.path.join(OUTPUT_DIR,"IndelsHistogram"))
os.makedirs(os.path.join(OUTPUT_DIR,"reports/AppendedReports"))
os.mkdir(os.path.join(GLM_DIR,"AppendGLM"))

subprocess.check_call("rm {LOG_DIR}/*".format(LOG_DIR),shell=True)
subprocess.check_call("rm {LOG_DIR}/MasterError.txt".format(LOG_DIR=LOG_DIR),shell=True)

subprocess.check_call("python {MACHETE}/writeStemIDFiles.py -o {ORIG_DIR} -f {OUTPUT_DIR}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR),shell=True)

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
popen = subprocess.Popen("wc -l {StemFile}".format(StemFile=StemFile),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
stdout,stderr = popen.communicate()
NUM_FILES = stdout
print(NUM_FILES)

## sorting reg files

#j1_id
processes = {}
for index in range(1,NUM_FILES + 1):
	cmd = "{MACHETE}/AlphabetizeENCODEreads.sh {ORIG_DIR}/reg/ {OUTPUT_DIR}} | awk '{print $4}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR)
	stdout = os.path.join(LOG_DIR,str(index) + "_out_1sortReg.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_1sortReg.txt")
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


print("sorting reg files")


# sorting genome files
#j2_id

print("sorting genome files")
processes = {}
for index in range(1,NUM_FILES + 1):
	cmd = "{MACHETE}/AlphabetizeENCODEreads.sh {ORIG_DIR}/genome {OUTPUT_DIR} | awk '{print $4}'".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR)
	stdout = os.path.join(LOG_DIR,str(index) + "_out_1sortGenome.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_1sortGenome.txt")
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
#
#
## finding mismatched paired end reads


#j3_id
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_2PEfinder.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_2PEfinder.txt")
	cmd = "{MACHETE}/PEfinder_genomeAndReg_ENCODE.sh {ORIG_DIR} {OUTPUT_DIR} {USERBPDIST} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR,USERBPDIST=USERBPDIST)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
	

echo "Outputting Mismatched paired ends: job ${j3_id}"

### window around position in genome is 10,000

### counting rates of mismatched PE

#
#j4_id
print("counting mismatch rates")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_3PEcounter.txt")
	stderr = os.path.joinn(LOG_DIR,str(index) + "_err_3PEcounter.txt")
	cmd = "{MACHETE}/DistantPE_Counter_genome_ENCODE.sh {OUTPUT_DIR} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# sort mismatched PE by chromosome

#j5_id
print("sorting mismatched PE files")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_4PEsort.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_4PEsort.txt")
	cmd = "{MACHETE}/SortPairedEnds.sh {OUTPUT_DIR} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#make Junction fasta file by extracting sequence info from pickles

print("make fusion fasta files")
processes = {}
for i in range(1,25):
	if i == 23:
		i = "X"
	elif i == 24:
		i = "Y"
	stdout = os.path.join(LOG_DIR,str(i) + "_out_5makefasta.txt")
	stderr = os.path.join(LOG_DIR,str(i) + "_err_5makefasta.txt")
	cmd = "{MACHETE}/makeJunctions.sh {EXONS} {OUTPUT_DIR} {i} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,EXONS=EXONS,OUTPUT_DIR=OUTPUT_DIR,i=i)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

##
##
##make single FJ fasta from all the fastas and then call bowtie indexer
##
#
#j6a_id
print("make FJ bowtie indices for each experiment")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_5FJIndexing.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_5FJIndexing.txt")
	cmd = "{MACHETE}/linkfastafiles.sh {OUTPUT_DIR} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
##
##
# make BadJunc directory --  bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align
#

genomeIndex = os.path.join(CIRCREF,"hg19_genome")
transcriptomeIndex = os.path.join(CIRCREF,"hg19_transcriptome")
regIndex = os.path.join(CIRCREF,"hg19_junctions_reg")
juncIndex = os.path.join(CIRCREF,"hg19_junctions_scrambled")


#j7_id=
print("Identify Bad FJ's")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_6BadJunc.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_6BadJunc.txt")
	cmd = "{MACHETE}/LenientBadFJ_SLURM.sh {OUTPUT_DIR} {REFGENOME} {MACHETE} {CIRCREF} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,REFGENOME=REFGENOME,CIRCREF=CIRCREF)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)



# align unaligned files to the FJ bowtie index
#
#j8_id
print("align unaligned reads to FJ index")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_7AlignFJ.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_7AlignFJ.txt")
	cmd = "{MACHETE}/AlignUnalignedtoFJ.sh {OUTPUT_DIR} {ORIG_DIR} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

##
#
#
###make FJ naive report
#
#j9_id
print("make naive rpt")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_8NaiveRpt.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_8NaiveRpt.txt")
	cmd = "{MACHETE}/FarJuncNaiveReport.sh {OUTPUT_DIR} {ORIG_DIR} {NUMBASESAROUNDJUNC} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,ORIG_DIR=ORIG_DIR,NUMBASESAROUNDJUNC=NUMBASESAROUNDJUNC)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

###
##
####### GLM #######
##Make Class input files for GLM
##
###
#j15a_id
print("make FJ class input files")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_15FJforGLM.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_15FJforGLM.txt")
	cmd = "{MACHETE}/parse_FJ_ID_for_GLM.sh {OUTPUT_DIR} {CIRCPIPE_DIR}/circReads/ | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


#
##
###
###### ESTIMATE LIGATION ARTIFACT
#
###
###make FarJunctions.fa => Indel.fa files
####
#j10_id
print("make indel files")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_10FJIndels.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_10FJIndels.txt")
	cmd = "{MACHETE}/MakeIndelFiles.sh {OUTPUT_DIR} {NumIndels} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,NumIndels=NumIndels)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

#
#
# make Bowtie Indices for Far Junc Indel files
#
print("index indels")
for i in range(1,NumIndels + 1):
	processes = {}
	for j in range(1,NUM_FILES + 1):
		stdout = os.path.join(LOG_DIR,"NumIndels{i}_{j}_out_11indexindels.txt".format(i=i,j=j))
		stderr = os.path.join(LOG_DIR,"NumIndels(i}_{j}_err_11indexindels.txt".format(i=i,j=j)
		cmd = "{MACHETE}/BowtieIndexFJIndels.sh {OUTPUT_DIR}/FarJuncIndels {i} {OUTPUT_DIR}/BowtieIndels {OUTPUT_DIR} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,i=i)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)


##Align FarJuncSecondary (unaligned to FJ index) to FJ indels


#
##
print("align to indels")
BOWTIEPARAMETERS="--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
for i in range(1,NumIndels + 1):
	processes = {}
	for j in range(1,NUM_FILES + 1):
		stdout = os.path.join(LOG_DIR,"NumIndels{i}_{j}_out_12alignindels.txt".format(i=i,j=j))
		stderr = os.path.join(LOG_DIR,"NumIndes{i}_{j}_err_12alignindels.txt".format(i=i,j=j))
		cmd = "{MACHETE}/BowtieAlignFJIndels.sh {OUTPUT_DIR} {BOWTIEPARAMETERS} {i} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,BOWTIEPARAMETERS=BOWTIEPARAMETERS,i=i)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)

#
## loop through AlignedIndels directory
## things that don't overlap indels junction by ${NumBPOverlapAtJunc} are removed.
## for every junction, a string is created.  For example, if 3 indel files exist, the string [0, 0, 2, 1, 5, 0, 1] represents 0*3 deletions, 0* 2 deletions, 2 * 1 deletion, 1 * no indels, 5 * 1 insertion, 0 * 2 insertions, 1 * 3 insertions.
## strings are output into IndelsHistogram folder
#
##
#j14_id
print("make indels histo")
processes = {}
for index in range(1,NumIndels + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_13filterIndels.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_13filterIndels.txt")
	cmd = "{MACHETE}/FindAlignmentArtifact_SLURM.sh {OUTPUT_DIR} {NumBPOverlapAtJunc} {NumIndels} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,NumIndels=NumIndels)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


###
####
##### REG INDELS ################
##
# Align unaligned files to the expanded reg junctions with indels
AlignedIndels = os.path.join(CIRCPIPE_DIR,"orig/RegIndelAlignments")
os.makedirs(AlignedIndels)

#j16_id
print("Aligning unaligned files to linear junc indels")
BOWTIEPARAMETERS = "--no-sq --no-unal --score-min L,0,-0.24 --rdg 50,50 --rfg 50,50"
for i in range(1,NumIndels + 1):
	processes = {}
	for j in range(1,NUM_FILES + 1):
	 	stdout = os.path.join(LOG_DIR,"NumIndels{i}_{j}_out_15AlignRegIndels.txt".format(i=i,j=j))
		stderr = os.path.join(LOG_DIR,"NumIndels{i}_{j}_err_15AlignRegIndels.txt".format(i=i,j=j))
		cmd = "{MACHETE}/AlignUnalignedtoRegIndel.sh {CIRCPIPE_DIR} {i} {OUTPUT_DIR} {BOWTIEPARAMETERS} {REG_INDEL_INDICES} | awk '{print $4}'".format(MACHETE=MACHETE,CIRCPIPE_DIR=CIRCPIPE_DIR,i=i,OUTPUT_DIR=OUTPUT_DIR,BOWTIEPARAMETERS=BOWTIEPARAMETERS,REG_INDEL_INDICES=REG_INDEL_INDICES)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
	checkProcesses(processes)

###
####
######## MAKE REG AND FJ INDELS CLASS OUTPUT FILES ###########
####
###
### reg indels class output file
#

#j18_id
print("Reg Indels Class Output")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_18RegIndelsClassOutput.txt")
	stdrr = os.path.join(LOG_DIR,str(index) + "_err_18RegIndelsClassOutput.txt")
	cmd = "{MACHETE}/RegIndelsClassID.sh {OUTPUT_DIR} {CIRCPIPE_DIR} {NumBPOverlapAtJunc} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,NumBPOverlapAtJunc=NumBPOverlapAtJunc)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# FJ indels class output file

#j19_id
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_19FJIndelsClassOutput.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_19FJIndelsClassOutput.txt")
	cmd = "{MACHETE}/FJIndelsClassID.sh {OUTPUT_DIR} {CIRCPIPE_DIR} {NumBPOverlapAtJunc} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,CIRCPIPE_DIR=CIRCPIPE_DIR,NumBPOverlapAtJunc=NumBPOverlapAtJunc)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
echo "FJ Indels Class Output: ${j19_id}"


###### RUN GLM ###########################
#
## Run GLM
##
#j15b_id
print("Run GLM")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_15GLM_r.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_15GLM_r.txt")
	cmd = "{MACHETE}/run_GLM.sh {CIRCPIPE_DIR} {OUTPUT_DIR} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,CIRCPIPE_DIR=CIRCPIPE_DIR,OUTPUT_DIR=OUTPUT_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)


## Append linear junctions GLM report with anomalies, indels
#
#j17_id
print("Appending linearJuncs GLM report")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_17AppendRegGLM.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_17AppendGLM.txt")
	cmd = "{MACHETE}/AddIndelstoGLM.sh {CIRCPIPE_DIR} {OUTPUT_DIR} {NumBPOverlapAtJunc} | awk '{print $4}'".format(MACHETE=MACHETE,CIRCPIPE_DIR=CIRCPIPE_DIR,OUTPUT_DIR=OUTPUT_DIR,NumBPOverlapAtJunc=NumBPOverlapAtJunc)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
##
####

### Append Naive report:  Add Indels, quality of junctions (good/bad), frequency of junction participation in linear junctions, and GLM report fields to the Naive report
###
#j15_id
print("append naive rpt")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_14AppendRpt.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_14AppendRpt.txt")
	cmd = "{MACHETE}/AppendNaiveRept.sh {OUTPUT_DIR} {GLM_DIR} {MACHETE} | awk '{print $4}'"
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
