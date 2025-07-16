# Load required library
library(data.table)

# ==== Config ====
dry_run <- TRUE  # Set to FALSE to actually copy/rename

# ==== Utility Functions with Dry Run Support ====

safe_copy <- function(from, to_dir, dry_run = FALSE) {
  if (length(from) == 0) {
    message("No files to copy.")
    return(invisible(NULL))
  }
  if (!dir.exists(to_dir) && !dry_run) dir.create(to_dir, recursive = TRUE)

  if (dry_run) {
    cat("DRY RUN: Would copy files to", to_dir, ":\n", paste(from, collapse = "\n"), "\n\n")
  } else {
    success <- file.copy(from = from, to = to_dir, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
    if (any(!success)) warning("Some files failed to copy: ", paste(from[!success], collapse = ", "))
  }
}

safe_rename <- function(old_names, new_names, dry_run = FALSE) {
  if (length(old_names) == 0) {
    message("No files to rename.")
    return(invisible(NULL))
  }
  if (dry_run) {
    cat("DRY RUN: Would rename files:\n")
    cat(paste(old_names, "→", new_names), sep = "\n")
    cat("\n")
  } else {
    success <- file.rename(from = old_names, to = new_names)
    if (any(!success)) warning("Some files failed to rename: ", paste(old_names[!success], collapse = ", "))
  }
}

# ==== 1. Copy .gz Files Matching Patterns ====

gz_files <- list.files('./', pattern = 'both.*\\.gz|host.*\\.gz')
safe_copy(gz_files, to_dir = '/PHShome/jn22/hERV/raw_data', dry_run = dry_run)

# ==== 2. Copy .fa Files ====

fasta_files <- list.files('./', pattern = '\\.fa$')
safe_copy(fasta_files, to_dir = '/PHShome/jn22/hERV/hERV_Work/', dry_run = dry_run)

# ==== 3. Rename .1. and .2. in .gz Files ====

gz_files <- list.files('./', pattern = '\\.gz$')
renamed1 <- sub('_1\\.', '_R1.', gz_files)
renamed1 <- sub('_2\\.', '_R2.', renamed1)
safe_rename(gz_files, renamed1, dry_run = dry_run)

# ==== 4. Copy First Mouse Fasta File ====

mouse_fa <- list.files('./', pattern = 'mouse.*\\.fa$')
if (length(mouse_fa) >= 1) {
  safe_copy(mouse_fa[1], to_dir = '../../scratch/erv/hERV_Work/', dry_run = dry_run)
} else {
  warning("No mouse fasta file found.")
}

# ==== 5. View First 100 Rows of RepeatMasker Files ====

if (!dry_run) {
  try({
    cb <- fread(nrows = 100, file = '/data/wanglab_mgberis/siyi2022summer_fastq2counts/hERV_Work/mouse_erv_combined_repeatMasker_on_m39.fa')
    # chr10 <- fread(nrows = 100, file = 'gtfs_for_fasta_needed_by_salmon/chr10.fa')
  }, silent = TRUE)
} else {
  cat("DRY RUN: Would read first 100 rows of RepeatMasker fasta files.\n\n")
}

# ==== 6. Final Rename .R1. → _R1_ and .R2. → _R2_ ====

gz_files <- list.files('./', pattern = '\\.gz$')
renamed2 <- sub('\\.R1\\.', '_R1_', gz_files)
renamed2 <- sub('\\.R2\\.', '_R2_', renamed2)
safe_rename(gz_files, renamed2, dry_run = dry_run)
