list.files('./',pattern = 'both.*.gz|host.*.gz')
file.copy(from=list.files('./',pattern = 'both.*.gz|host.*.gz'),
          to='/PHShome/jn22/hERV/raw_data', 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)

list.files('./',pattern = '.fa$')
file.copy(from=list.files('./',pattern = '.fa$'),
          to='/PHShome/jn22/hERV/hERV_Work/', 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)

list.files('./',pattern = '.gz$')->curren_name
new_name=sub(pattern = '\\.1\\.',replacement = '\\.R1\\.',x = curren_name)
new_name=sub(pattern = '\\.2\\.',replacement = '\\.R2\\.',x = new_name)

file.rename(curren_name, 
            new_name)

list.files('./',pattern = 'mouse.*fa$')->fl
file.copy(from=fl[1],
          to='../../scratch/erv/hERV_Work/', 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)

# there were repeat IDs in the combined fasta
library(data.table)
setwd("/PHShome/jn22/siyi_summer2023/erv/00repeatmasker")
cb<-fread(nrows = 100,file = '/PHShome/jn22/hERV/hERV_Work/mouse_erv_combined_repeatMasker_on_m39.fa')
chr10<-fread(nrows = 100,
                   file = '/PHShome/jn22/siyi_summer2023/erv/00repeatmasker/gtfs_for_fasta_needed_by_salmon/chr10.fa')


list.files('./',pattern = '.gz$')->curren_name
new_name=sub(pattern = '\\.R1\\.',replacement = '\\_R1\\.',x = curren_name)
new_name=sub(pattern = '\\.R2\\.',replacement = '\\_R2\\.',x = new_name)

file.rename(curren_name, 
            new_name)
