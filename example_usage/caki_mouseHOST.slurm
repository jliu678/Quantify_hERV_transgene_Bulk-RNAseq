#!/bin/bash

### the script succeeds from test_single.lsf

#SBATCH --partition=normal
#SBATCH --job-name=ut
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=72G
#SBATCH --output=/data/wanglab_mgberis/__output.%j_%x
#SBATCH --error=/data/wanglab_mgberis/__error.%j_%x
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jliu78@mgh.harvard.edu

# specify the path to the folder contains both herv_project/ and hERV_Work/ and raw_data/

working_dir="/data/wanglab_mgberis/siyi2022summer_fastq2counts"

cd "${working_dir}/herv_project"

if [ ! -d "logs" ]; then
    mkdir logs
fi

echo "Job started at: $(date)"
echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID}"

### for mgh cluster, keep -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1

bash main.sh -PLATFORM "cluster-mgh" -ANALYSIS_STEP "all" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 -REF_GENOME "https://a/GRCm39.primary_assembly.genome.fa" \
-TRANSCRIPTS "https://a/mouse_erv_combined_repeatMasker_on_m39.fa" \
-TRANSCRIPTS "transgene.fa" \
-REF_TRANSCRIPT "https://a/mouse_gencode.vM33.transcripts.fa.gz"

### for local desktop like the imagedeskto, either not supply '-PLATFORM' parameter or supply anything else, for example
### bash main.sh -PLATFORM "anything_different_with_cluster-mgh"

################## OUTPUT FILE ######################################
### all output files will be stored in ../hervwork/results

echo "Job completed at: $(date)"