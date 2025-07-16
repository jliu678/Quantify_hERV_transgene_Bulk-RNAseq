# 🧬Quantify hERV and transgenes from Bulk-RNAseq

> 💡 **Tip: [Please find **👉 MY BLOG** for an introduction to the project, along with the detailed mathematical and biological reasoning behind the code in this repository.](https://myhugoblog)**


## 📑 Table of Contents
- [📘 Introduction & Reasoning](#-introduction--reasoning-complete-version-is-here)
- [💡 Usage Example](#-usage-example)
  - [💻 LSF Scheduler on High-performance Cluster](#-lsf-scheduler-on-high-performance-cluster)
    - [1. Configure LSF Parameters and Working Directory](#1-configure-lsf-parameters-and-working-directory)
    - [2. Wrangle File Input](#2-wrangle-file-input)
    - [3. Other Parameters](#3-other-parameters)
    - [4. Output File](#4-output-file)
  - [🧾 Other Usage Examples](#-other-usage-examples)

## 📘Introduction & Reasoning (complete version is [here](myhugoblod))

Briefly, we use `Salmon` to quantify the hERV and transgenes from bulk-RNAseq because of the (Bayesian) EM algorithm implenmented in `Salmon`. EM algorithm is advantageous in dealing with  multimaping, which is commonly seen for hERV and transgene quantification, because:
- mathematically, EM is well suitable for estimating parameters in Gaussian Mixture Models, which is similar to the model describing multimapping
- biologically, EM reasonably takes into consideration the mapping probability and mapping bias

The source codes here can be deploied in Cloud Cluster Computational platform and local desktops, and allows for: 
- batch computaion with tunable batch size
- automatic environment setting up
- file preprocessing
- QC
- index buidling, alignment and feature counting by either `subread` or `salmon`

The github repo also contains [usage examples using `SLURM` or `LSF` job scheduler](https://github.com/jliu678/herv_project_siyi/tree/main/example_usage) and [handy utility tools](https://github.com/jliu678/herv_project_siyi/tree/main/utils).

## 💡usage example

### 💻 LSF scheduler on High-performance Cluster

#### 1. configure LSF parameters and working directory

```bash
if [ ! -d "logs" ]; then
	mkdir logs
fi 
#BSUB -e logs/%J-errors.err
#BSUB -o logs/%J-outputs.out

#BSUB -q normal
#BSUB -n 16
#BSUB -R 'rusage[mem=32000]'

cd ~/hERV/herv_project
```

#### 2. Wrangle file input

- Stores in `../raw_data` all fastq files named as 'fq' or better 'fq.gz', and accordingly specify `-SOURCE raw_data` of main.sh

- DO NOT use `random_name.r1.fq.gz` as raw fastq file name, but instead use `random_name_r1.fq.gz`, because he postfix of the raw fastq file names `$name` is obtained by `${name#*.}`. Please see below cases:

   ```bash
   (base)$ name=random_name.fq.gz; echo ${name#*.}
     fq.gz
   (base)$ name=random_nam.e.fq.gz; echo ${name#*.}
     e.fq.gz
   ```
- For subread, store in `../hERV_Work` all the annotation gtf/gff3 files and fasta files 

- For salmon, if you need generate index, please specify below parameters:

    - `-REF_GENOME`: to specify the genome fasta used to generate decoys
    - `-REF_TRANSCRIPT` and `TRANSCRIPTS` are transcripts to count, which are needed for indexing and counting
    
- For salmon, if you need generate index, first option is to supply parameter `-REF_GENOME`, `-REF_TRANSCRIPT` and `TRANSCRIPTS` of "main.sh" with specific link, for example:
   ```bash
   bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
   -REF_GENOME "https://.*/weird.genome.fa.gz" \
   -REF_TRANSCRIPT "https://.*/transcripts_of_widely_Used_Normal_OR_Commen_Genes.fa.gz" \
   -TRANSCRIPTS "https://.*/transcripts_of_erv_or_other_special_features.fa.gz" \
   -TRANSCRIPTS "https://.*/if_you_have_another_special_features.fa.gz"
  ```

- For salmon, if you need generate index, 2nd option is to store the already downloaded "fa.gz" or "fa" files in .`./hERV_Work/` to avoid downloading, their file names need follow below convention: 
 
   (A) genome fasta can be named as 'GRCh38.p13.genome.fa' which matchs with `REF_GENOME_LOC='https://.*/GRCh38.p13.genome.fa.gz'` in main.sh. Internal functions called will get a string by removing everything before the last "/" and the last "/" itself and remove ".gz", and if there is no file titled same as the string, it will download from the http. 

   (B) genome fasta can be named as something else but make sure you specify parameter `-REF_GENOME` of "bash main.sh" with corresponding valuechange (the argument is passed to the `REF_GENOME_LOC` variable in main.sh), for example if the fasta file is `weird.genome.fa`:

   ```bash
   bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
   -REF_GENOME "weird.genome.fa"
   ```                                     

   (C) similarly, reference genome fasta can be named as 'gencode.v41.transcripts.fa.gz' which matchs with `REF_TRANSCRIPT_LOC='https://.*/gencode.v41.transcripts.fa.gz'` in main.sh;  or name something else but make sure you change the value of the `REF_TRANSCRIPT_LOC` variable in main.sh by supplying parameter `-REF_TRANSCRIPT` of "bash main.sh" with corresponding value, for example to match you fasta file you name as `weird.transcripts.fa`:

   ```bash
   bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
   -REF_TRANSCRIPT "/weird.transcripts.fa.gz"
   ```
   (D) likewise, other fasta can be named as `package-entities-erv.fa` which matchs with `TRANSCRIPT_LOCS='https://.*/package-entities-erv.fa.gz'` in main.sh; or named as something else but make sure you change the value of the `TRANSCRIPT_LOCS` variable in main.sh by supplying parameter `-TRANSCRIPTS` of "bash main.sh" with corresponding value, for example to match you special fq.gz file you name as `ugly.transcripts.fa.gz`:
   
   ```bash
   bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
   -TRANSCRIPTS "https://.*/ugly.transcripts.fa.gz" 
   ```
- For salmon, if you have index already, you can bypass re-generating index by storing the index folder in `../hERV_Work/salmon/`, and rename the index folder as the output of below scripts:

    ```bash
    
    #-TRANSCRIPTS) TRANSCRIPT_LOCS+=("$2")
    
    TRANSCRIPTS=$(basename $(basename $TRANSCRIPT_LOCS .gz) .bz2)
    
    echo $(IFS=-; echo "${TRANSCRIPTS[*]%.*}")
    
    ```  
    where `TRANSCRIPT_LOCS` is an array to contain all the values you will assign to parameter `-TRANSCRIPTS`. For example, 
    ```bash
    TRANSCRIPTS=(
    "https://.*/weird.transcripts.fa.gz"
    "https://.*/2nd.transcripts.fa.gz"
    "https://.*/3rd.fa.gz"
    )
    ```
    for 
    ```bash
    bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
   -TRANSCRIPTS "https://.*/weird.transcripts.fa.gz"
   -TRANSCRIPTS "https://.*/2nd.transcripts.fa.gz"
   -TRANSCRIPTS "https://.*/3rd.fa.gz"
    ```

- If you want to input multiple fasta that are not reference genome, you have to supply multiple times of parameter `-TRANSCRIPTS` to "bash main.sh", for example you can store already downloaded `weird.transcripts.fa`, `2nd.transcripts.fa.gz` and `3rd.fa.gz` in `../hERV_Work` and speficy:

   ``` bash 
   bash main.sh -PLATFORM "cluster-mgh" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 \
   -TRANSCRIPTS "https://.*/weird.transcripts.fa.gz"
   -TRANSCRIPTS "https://.*/2nd.transcripts.fa.gz"
   -TRANSCRIPTS "https://.*/3rd.fa.gz"
   ```

- The script will identify files as paired when there are two files named in the way that there is file2 whose name is either `sub(name_of_file1, pattern="r1", replacement="r2")` or `sub(name_of_file1, pattern="R1", replacement="R2")`. [`sub` here is base R function](https://stat.ethz.ch/R-manual/R-devel/library/base/html/grep.html). See `group_fastq()` and `check_name()` in "analysis.sh" for detail.


#### 3. other parameters
- parameter `EXIT_ON_SINGLE` of main.sh has default value of `false`, thus single-end sequenced fq files can be analysis using the same bash command

- for mgh cluster, please keep `-BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1`, for example:

   ```bash
   bash main.sh -PLATFORM "cluster-mgh" -ANALYSIS_STEP "all" -SOURCE "raw_data" -BATCH_SIZE 1 -CLEAR_TMP -MAX_PARALLEL 1 
   ```
- for local desktop like the imagedeskto, either not supply `-PLATFORM` parameter or supply anything else, for example

   ```bash
   bash main.sh -PLATFORM "anything_different_with_cluster-mgh"
   ```
#### 4.  output file
 - all output files will be stored in ../hervwork/results

### 🧾 other usage examples
Please see [usage examples using `SLURM` or `LSF` job scheduler](https://github.com/jliu678/herv_project_siyi/tree/main/example_usage)
