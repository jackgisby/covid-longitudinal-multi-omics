##Section 1: Header                                                    #
#########################################################################

#!/bin/bash
#PBS -N ArrayJob
#PBS -l select=2:ncpus=16:mem=62gb
#PBS -l walltime=08:00:00
#PBS -m abe
#PBS -o jobout
#PBS -e joberr

#PBS -J 1-378
#

#########################################################################
## Section 2: Transferring Data from current directory to temp folder in $EPHEMERAL/job_array     #
#########################################################################

# Name the Project directory

PROJECT_SUBDIR="ArrayJob"
JOB_ID="HTseq-Gencode"

# Specify the filenames
#INPUTFILES=(C139-3001_R1_001.fastq.gz  C168-2801_R2_001.fastq.gz  C23-2204_R1_001.fastq.gz)

INPUTFILES_DIRECTORY="$RDS_PROJECT/peters-lab-sle-multiomics/live/temp_sle_covid19_transcriptomics/star_salmon/"

INPUTFILES=($(ls $INPUTFILES_DIRECTORY*.bam))

# Pull data file name from list defined above 
#
INPUTFILEPATH="${INPUTFILES[$PBS_ARRAY_INDEX-1]}"

INPUTFILENAME=$(basename "${INPUTFILEPATH%.*}")

# Specify working directory under the /temp
WORKING_DIR="$EPHEMERAL/job_array/$PROJECT_SUBDIR-$JOB_ID"

# If working directory does not exist, create it
# The -p means "create parent directories as needed"
if [ ! -d "$WORKING_DIR" ]; then
    mkdir -p $WORKING_DIR
fi

# Specify destination directory (this will be subdirectory of your user directory in the archive)
DESTINATION_DIR="/rds/general/project/peters-lab-sle-multiomics/live/temp_sle_covid19_transcriptomics/$JOB_ID"

# If destination directory does not exist, create it
 #The -p in mkdir means "create parent directories as needed"
if [ ! -d "$DESTINATION_DIR" ]; then
    mkdir -p $DESTINATION_DIR
fi
#
# Copy the input data to the working directory
#cp $INPUTFILENAME $WORKING_DIR/



#########################################################################
## Section 3:Executing the program                                      #
#########################################################################

# navigate to the working directory

cd $WORKING_DIR

###test
#

#ls "$INPUTFILEPATH" > "$INPUTFILENAME".txt

# Run the program

module load samtools

module load anaconda3/personal

source activate HTseq

GTFFILE="/rds/general/project/peters-lab-sle-multiomics/live/temp_sle_covid19_transcriptomics/genome/gencode.v38.annotation.gtf"

python -m HTSeq.scripts.count --mode=intersection-strict --format=bam --stranded=reverse --type=gene $INPUTFILEPATH $GTFFILE > $INPUTFILENAME.count


#########################################################################
## Section 4: Copy the result to /archive                               #
#########################################################################

scp * $DESTINATION_DIR

# clear the working directory
##rm -rf $WORKING_DIR
