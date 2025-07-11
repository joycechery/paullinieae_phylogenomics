#!/bin/bash
#SBATCH --job-name=00.trimmomatic_batch3.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=32GB
#SBATCH --time=40:00:00
#SBATCH --output=job.%A_%a.out
#SBATCH --error=job.%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

# Load Trimmomatic module if required by your system
module load trimmomatic/0.39
module load openjdk/21.0.3

# Define paths to Trimmomatic, adapter file, and input/output directories
TRIMMOMATIC_JAR="/share/apps/trimmomatic/0.39/trimmomatic-0.39.jar"
ADAPTERS="/share/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa"
INPUT_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/raw_data/batch3"
OUTPUT_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/00.trimmomatic/clean_data/batch3"
UNPAIRED_DIR="${OUTPUT_DIR}/unpaired"

# Ensure output and unpaired directories exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$UNPAIRED_DIR"

# Path to the namelist file
NAMELIST="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/scripts/namelist/namelist_batch3.txt"

# Read and process each name in the namelist
while IFS= read -r BASE || [[ -n $BASE ]]; do
  # Define input files
  FORWARD_READ="${INPUT_DIR}/${BASE}_R1.fastq.gz"
  REVERSE_READ="${INPUT_DIR}/${BASE}_R2.fastq.gz"

  # Define output trimmed files
  FORWARD_PAIRED="${OUTPUT_DIR}/${BASE}_R1_paired.fastq.gz"
  FORWARD_UNPAIRED="${OUTPUT_DIR}/${BASE}_R1_unpaired.fastq.gz"
  REVERSE_PAIRED="${OUTPUT_DIR}/${BASE}_R2_paired.fastq.gz"
  REVERSE_UNPAIRED="${OUTPUT_DIR}/${BASE}_R2_unpaired.fastq.gz"

  # Check if the output files already exist
  if [[ -f $FORWARD_PAIRED && -f $REVERSE_PAIRED ]]; then
    echo "Output files for $BASE already exist. Skipping..."
    continue
  fi

  # Run Trimmomatic
  echo "Processing $BASE..."
  java -jar "$TRIMMOMATIC_JAR" PE \
    -phred33 \
    "$FORWARD_READ" "$REVERSE_READ" \
    "$FORWARD_PAIRED" "$FORWARD_UNPAIRED" \
    "$REVERSE_PAIRED" "$REVERSE_UNPAIRED" \
    ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36

  echo "$BASE processing completed."
done < "$NAMELIST"

# Move all _unpaired.fastq.gz files to the unpaired folder
echo "Organizing unpaired files..."
find "$OUTPUT_DIR" -name "*_unpaired.fastq.gz" -exec mv {} "$UNPAIRED_DIR" \;

# Delete the unpaired folder and its contents
echo "Deleting unpaired files..."
rm -rf "$UNPAIRED_DIR"

# Rename all _paired.fastq.gz files to _clean.fastq.gz
echo "Renaming paired files..."
find "$OUTPUT_DIR" -name "*_paired.fastq.gz" | while read -r FILE; do
  NEW_NAME="${FILE/_paired.fastq.gz/_clean.fastq.gz}"
  mv "$FILE" "$NEW_NAME"
done

echo "All files processed and organized."
