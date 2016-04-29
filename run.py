import os
import sys
import subprocess
import glob
from argparse import ArgumentParser

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


ORIG_DIR = os.path.join(CIRCPIPE_DIR,"orig")
UNALIGNEDDIR = os.path.join(ORIG_DIR),"unaligned")
GLM_DIR = os.path.join(CIRCPIPE_DIR,"circReads/glmReports")
DistantPEDir = os.path.join(OUTPUT_DIR,"DistantPEFiles")
os.mkdir(DistantPEDir)
FASTADIR = os.path.join(OUTPUT_DIR,"fasta")
os.mkdir(FASTADIR)
BOWTIE_DIR = os.path.join(OUTPUT_DIR,"BowtieIndex")
os.mkdir(BOWTIE_DIR)
FARJUNCDIR = os.path.join(OUTPUT_DIR,"FarJunctionAlignments")
os.mkdir(FARJUNCDIR)
SECONDFARJUNCDIR = os.path.join(OUTPUT_DIR,"FarJuncSecondary")
os.mkdir(SECONDFARJUNCDIR)
BadFJDir = os.path.join(OUTPUT_DIR,"BadFJ")
os.mkdir(BadFJDir)
StemFile = os.path.join(OUTPUT_DIR,"StemList.txt")

os.mkdir(os.path.join(OUTPUT_DIR,"reports"))
LOG_DIR = os.path.join(OUTPUT_DIR,"err_and_out")
os.mkdir(LOG_DIR)
os.mkdir(os.path.join(OUTPUT_DIR,"BowtieIndels"))
os.mkdir(os.path.join(OUTPUT_DIR,"FarJuncIndels"))
os.mkdir(os.path.join(SECONDFARJUNCDIR,"AlignedIndels"))
os.mkdir(os.path.join(OUTPUT_DIR,"IndelsHistogram"))
os.mkdir(os.path.join(OUTPUT_DIR,"reports/AppendedReports"))
os.mkdir(os.path.join(GLM_DIR,"AppendGLM"))
os.mkdir(os.path.join(OUTPUT_DIR,"GLM_classInput"))

subprocess.check_call("rm {LOG_DIR}/*".format(LOG_DIR),shell=True)
subprocess.check_call("rm {LOG_DIR}/MasterError.txt".format(LOG_DIR=LOG_DIR),shell=True)

## This python script detects all the unique names for all pairs of files within a directory, eg. SRR12345, SRR123456, etc into a file called ${StemFile}
if os.path.isfile(StemFile):
	print("using existing StemList.txt")
else:
	print("generating StemList.txt from KNIFE output directory filenames")
	subprocess.check_call("python {MACHETE}/writeStemIDFiles.py -o {ORIG_DIR} -f {OUTPUT_DIR}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR),shell=True)

# counting # of times to go through the "PE matching" step - is the number of paired genome files ending in .sam /2
NUM_FILES = len(open(StemFile,"r").readlines())
NUM_FILES = stdout
print(NUM_FILES)

## if the program has been run before, there will be "sorted" reg and genome files.
## these are removed if they are present.
## All files from the original KNIFE alignments are sorted into alphabetical order because it is faster for python to identify read partners in two alphabetically sorted files than it is to find read pairs in two huge files where read IDs are scrambled.


## the shell AlphabetizeKNIFEreads.sh takes directories reg and genome, where we plan to search for mismatched paired ends, and sorts them alphabetically using the linux sort function

## sorting reg files
#j1_id
processes = {}
for index in range(1,NUM_FILES + 1):
	cmd = "{MACHETE}/AlphabetizeENCODEreads.sh {ORIG_DIR}/reg {OUTPUT_DIR} | awk '{print $4}".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR)
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
## The shell PEfinder.sh takes the KNIFE alignment directory (OrigDir -1 ), the output directory (FJDir -2 ), the distance beyond which the user would consider alignments to be "discordant" (BP_distance -3 ), and the MACHETE installation directory (4) and calls a python script PEfinder_genomeAndReg_ENCODE.py.
## The python script identifies paired R1 and R2 from the genome and reg alignments and if the alignment location is > user defined bases apart then records them in an output file within the output directory: FarJunctionDirectory/DistantPEFiles/<STEM>_distant_pairs.txt
## If, for example, a read pair was found to be discordant, and R1= chrA:X, R2=chrB:Y, then the distant_pairs.txt file would contain the readID and chrA:M-N, chrB:P-Q where M-N is a window of 10,000 bases on each side of X and P-Q is a window of 10,000 bases on each side of Y.
## The window of 10,000 bases can be set by the user in the shell PEfinder.sh
#j3_id
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_2PEfinder.txt")
	stderr = os.path.join(LOG_DIR,str(index) + "_err_2PEfinder.txt")
	cmd = "{MACHETE}/PEfinder.sh {ORIG_DIR} {OUTPUT_DIR} {USERBPDIST} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,ORIG_DIR=ORIG_DIR,OUTPUT_DIR=OUTPUT_DIR,USERBPDIST=USERBPDIST)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)
	

print("Outputting Mismatched paired ends: job ${j3_id}")
## Because there are lot of repeat locations in the _distant_pairs.txt file generated above, the DistantPE_Counter.py script is called by the shell of the same name to 1) eliminate duplicate locations.  Another early problem was that the fasta generated later was too enormous to run in a timely fashion so the distant_pairs.txt file is split by this script into 24 smaller files based on the chromosome # of the upstream partner.
## The shell DistantPE_Counter_genome_ENCODE.sh takes in the FarJunction output directory and MACHETE installation directory and outputs <FJDir>/DistantPEFiles/<STEM>/chr1,2,3,4,...,X,Y_Distant_PE_frequency.txt
#  The chrA_Distant_PE_frequency.txt files contain three columns: chrA:M-N, chrB:P-Q, and R, where R is the number of times that these two exact windows were matched together.  R could be used to cull the fasta file if it gets too large, but at this point we are still looking for junctions between exons if only one read pair aligned discordantly.
#
#j4_id
print("counting mismatch rates")
processes = {}
for index in range(1,NUM_FILES + 1):
	stdout = os.path.join(LOG_DIR,str(index) + "_out_3PEcounter.txt")
	stderr = os.path.joinn(LOG_DIR,str(index) + "_err_3PEcounter.txt")
	cmd = "{MACHETE}/DistantPE_Counter.sh {OUTPUT_DIR} {MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	processes[popen] = {"stdout":stdout,"stderr":stderr,"cmd":cmd}
checkProcesses(processes)

# sort mismatched PE by chromosome
## This is a simple shell script SortPairedEnds.sh to sort the chrA_Distant_PE_frequency.txt files into alphabetical order.  It takes FJDir/DistantPEFiles/chrA_Distant_PE_frequency.txt and outputs to same directory, sorted_chrA_Distant_PE_frequency.txt using the linux "sort" command.
## The reason for sorting is to increase the speed of the next step.

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
## Pickle files are binary files used for storage of large amounts of information in Python.  Here, GTF files have been repackaged into pickle files and store sequence, exon name, and exon location information.  Accessing pickles can be time consuming so target locations from the discordant reads have been broken down by chromosome and alphabetized.
## Loop goes through integers 1-24 where that corresponds to chromosome # (23 = X, 24 = Y). For each of those chromosomes, the shell makeJunctions.sh calls makeJunctions.py
## MakeJunctions.py takes in FJDir/DistantPEFiles/sorted__chrA_Distant_PE_frequency.txt, and outputs FJDir/fasta/<STEM>/<STEM>_chrAFarJunctions.fa.  It reads in the discordant windows, and searches the pickle file for all exons names/exon locations/ exon sequences on either the sense or antisense strands within the discordant windows.  Then it makes all pairs of possible exon junctions. All sequences are 300 base pairs long - 150 bases on each side of the breakpoint.  For exons that are fewer than 150 bases, the remainder of the 150 bases is padded with N's.  All pairs include fusions between two exons that are +/+, +/-, -/+, and -/-.  In the case that a sense and antisense exon are fused, then the sequence listed is the exact sequence that would occur if the fusion was read from the 5'->3' direction.  If the generated sequence was BLATted, the correct "strands" would appear in BLAT.  Similarly, for (-)/(-) exon pairs, if the generated sequence was BLATted, the exons would appear on the (-) strand in BLAT.
## THIS IS DIFFERENT IN KNIFE.  In KNIFE, if a (-)/(-) pair was BLATted, it would look as if it were +/+ because the KNIFE reverse complements (-) sequences.  Additionally in KNIFE, there is no way to detect inversions because -/+ and +/- fasta sequences are not generated.
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


##make single FJ fasta from all the fastas and then call bowtie indexer
##
## For each experiment, fasta files are generated for each chromosome separately as above.  The Bowtie2 call converts these into binary index files so the chromosome specific files must be concatenated into a single fasta file before generation of this index.
## The script linkfastafiles.sh uses linux to concatenate the <FJDir>/fasta/<STEM>/<STEM>_chr1,2,3,...,X,Y_FarJunctions.fa into a single large fasta <FJDir>/fasta/<STEM>_FarJunctions.fa.
## The second step of the linkfastafiles.sh calls Bowtie to build the Far Junctions bowtie index named <FJDir>/BowtieIndex/<STEM>/<STEM>_FJ_Index
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

## If there is homology between a FarJunctions fasta sequence and the genome or transcriptome or a linear junction or circular junction, then the fusion read is less likely.  Alignments of the FarJunctions fasta sequences to the KNIFE reference indices, genome, transcriptome, linear junctions (reg), and scrambled junctions (junc) are created with two different bowtie parameters.  Bad juncs will align to genome/transcriptome/junc/reg but good juncs will not align. These are just aligning the FJ Fasta to the bad juncs with various alignment parameters. Any junctions aligning to here will eventually be tagged as "BadFJ=1" in the final reports whereas if junctions don't align, they will receive a "BadFJ=0" in the final reports.

genomeIndex = os.path.join(CIRCREF,"hg19_genome")
transcriptomeIndex = os.path.join(CIRCREF,"hg19_transcriptome")
regIndex = os.path.join(CIRCREF,"hg19_junctions_reg")
juncIndex = os.path.join(CIRCREF,"hg19_junctions_scrambled")

# for BadFJ we Align FarJunc fasta file to the above indices with the following bowtie parameters:
# A minimum alignment score corresponding to 4 mismatches per 100 base pairs, no N ceiling, and a prohibitive read gap penalty that disallows any read gaps in the fasta sequence or the reference index.  Alignments are found in <FJDir>/BadFJ/<STEM>/<STEM>_BadFJto<ReferenceIndex>.sam.
BOWTIEPARAM="-f --no-sq --no-unal --score-min L,0,-0.24 --n-ceil L,0,100 -p 4 --np 0 --rdg 50,50 --rfg 50,50"

for i in range(1,NUM_FILES + 1):
	stemCmd = "awk 'FNR == '{c}' {{print $1}}' {StemFile}".format(i=i,StemFile=StemFile)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	stdout,stderr = popen.communicate()	
	retcode = popen.returncode
	if retcode:
		raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}."format(cmd=stemCmd,retcode=retcode,stdout=stdout,stderr=stderr))
	STEM = stdout
	FarJuncFasta = glob.glob(os.path.join(FASTADIR,"STEM","*FarJunctions.fa"))
	BadFJStemDir =os.path.join(BadFJDir,STEM)
	os.mkdir(BadFJStemDir)
	BadFJver2Dir=os.path.join(OUTPUT_DIR,"BadFJ_ver2",STEM)
	os.makedirs(BadFJver2Dir)
	r1file=os.path.join(BadFJver2Dir,STEM + "_FarJunctions_R1.fa")
	r2file=os.path.join(BadFJver2Dir,STEM + "_FarJunctions_R2.fa")
	#Prior to alignment with the reference indices a python script SplitFastaforBadFJ.py called by the shell LenientBadFJ_SLURM
	# is used to 1) remove all N's from the fasta sequences and 2) split the fasta sequence into a "read1" and "read2" file 
	# -- <FJdir>/BadFJ_ver2/<Stem>/<Stem>_FarJunctions_R1/2.fa.  The read 1s are the first 40 non-N bases and the read 2's are the 
	# last 40 non-N reads from the sequence.
	
	#j7_id=
	print("Identify Bad FJ's")
	processes = {}
	stdout = os.path.join(LOG_DIR,str(STEM) + "_out_6BadJunc.txt")
	stderr = os.path.join(LOG_DIR,str(STEM) + "_err_6BadJunc.txt")
	cmd = "{MACHETE}/LenientBadFJ_SLURM.sh ${FarJuncFasta} ${BadFJver2Dir} ${OUTPUT_DIR} ${MACHETE} | awk '{print $4}'".format(MACHETE=MACHETE,OUTPUT_DIR=OUTPUT_DIR,REFGENOME=REFGENOME,CIRCREF=CIRCREF)
	popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
	stdout,stderr = popen.communicate()

	BadFJtoGenomeFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoGenome.sam"):
	if os.path.exists(BadFJtoGenomeFile):
		print("{BadFJtoGenomeFile} exists. To realign, please manually delete this file first".format(BadFJtoGenomeFile=BadFJtoGenomeFile))
	else:
		stdout = os.path.join(BadFJStemDir,"out.txt")
		stderr = os.path.join(BadFJStemDir,"err.txt")
		cmd = "{MACHETE}/BowtieAligner.batch.sh {BOWTIEPARAM} {genomeIndex} {SPORKFasta} {BadFJtoGenomeFile} | awk '{print $4}'".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,SPORKFasta=SPORKFasta,BadFJtoGenomeFile=BadFJtoGenomeFile)
		print("BadFJ to genome: ${BadFJj1_id}")
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout,stderr = popen.communicate()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}."format(cmd=stemCmd,retcode=retcode,stdout=stdout,stderr=stderr))

	BadFJtotranscriptomeFile = os.path.join(BadFJStemDir,STEM + "__BadFJtotranscriptome.sam")
	if os.path.exists(BadFJtotranscriptomeFile):
		print("{BadFJtotranscriptomeFile} exists.  To realign, please manually delete this file first".format(BadFJtotranscriptomeFile=BadFJtotranscriptomeFile))
	else:
		stdout = os.path.join(BadFJStemDir,"out.txt")
		stderr = os.path.join(BadFJStemDir,"err.txt")
		cmd = "{MACHETE}/BowtieAligner.batch.sh {BOWTIEPARAM} {transcriptomeIndex} {SPORKFasta} {BadFJtotranscriptomeFile} | awk '{print $4}'".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,transcriptomeIndex=transcriptomeIndex,SPORKFasta=SPORKFasta,BadFJtotranscriptomeFile=BadFJtotranscriptomeFile)
		print("BadFJ to transcriptome")
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout,stderr = popen.communicate()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}."format(cmd=stemCmd,retcode=retcode,stdout=stdout,stderr=stderr))

	BadFJtoRegFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoReg.sam")
	if os.path.exists(BadFJtoRegFile):
		print("{BadFJtoRegFile} exists. To realign, please manually delete this file first.".format(BadFJtoRegFile=BadFJtoRegFile))
	else:
		stdout = os.path.join(BadFJStemDir,"out.txt") 
		stderr = os.path.join(BadFJStemDir,"err.txt")
		cmd = "{MACHETE}/BowtieAligner.batch.sh {BOWTIEPARAM} {regIndex} {SPORKFasta} {BadFJtoRegFile} | awk '{print $4}'".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,regIndex=regIndex,SPORKFasta=SPORKFasta,BadFJtoRegFile=BadFJtoRegFile)
		print("BadFJ to reg")
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout,stderr = popen.communicate()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}."format(cmd=stemCmd,retcode=retcode,stdout=stdout,stderr=stderr))


	
	BadFJtoJuncFile = os.path.join(BadFJStemDir,STEM + "_BadFJtoJunc.sam")
	if os.path.exists(BadFJtoJunc):
		print("{BadFJtoJuncFile} exists. To realign, please manually delete this file first".format(BadFJtoJuncFile=BadFJtoJuncFile))
	else
		stdout = os.path.join(BadFJStemDir,"out.txt")
		stderr = os.path.join(BadFJStemDir,"err.txt")
		print("BadFJ to junc: ")
		cmd = "{MACHETE}/BowtieAligner.batch.sh {BOWTIEPARAM} {juncIndex} {SPORKFasta} {BadFJtoJuncFile} | awk '{print $4}'".format(MACHETE=MACHETE,BOWTIEPARAM=BOWTIEPARAM,juncIndex=juncIndex,SPORKFasta=SPORKFasta,BadFJtoJuncFile=BadFJtoJuncFile)
		popen = subprocess.Popen(cmd,stdout=stdout,stderr=stderr,shell=True)
		stdout,stderr = popen.communicate()
		retcode = popen.returncode
		if retcode:
			raise Exception("Command {cmd} failed with return code {retcode}. stdout is {stdout} and stderr is {stderr}."format(cmd=stemCmd,retcode=retcode,stdout=stdout,stderr=stderr))










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
