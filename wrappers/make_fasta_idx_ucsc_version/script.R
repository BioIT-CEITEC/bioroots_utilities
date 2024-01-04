#########################################
# wrapper for rule: make_fasta_idx_ucsc_version
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[["run"]]
sink(logfile, append = T, type = "output")
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: make_fasta_idx_ucsc_version \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(GenomeInfoDb)
library(data.table)

### IN PROGRESS
## prepared manually with the use of GenomeInfoDb::genomeStyles
#organism_list <- c("Arabidopsis_thaliana","Caenorhabditis_elegans","Canis_familiaris","Cyanidioschyzon_merolae","Drosophila_melanogaster","Homo_sapiens","Mus_musculus","Oryza_sativa","Populus_trichocarpa","Rattus_norvegicus","Saccharomyces_cerevisiae","Zea_mays")
#names(organism_list) <- c("arabidopsis","name1","canis_familiaris","name2","drosophila_melanogaster")

input_file <- snakemake@input[["idx"]]
output_file <- snakemake@output[["ucsc"]]
dest_DB <- snakemake@params[["dest"]]
organism <- snakemake@params[["organism"]]

if(organism %in% c("arabidopsis","Arabidopsis thaliana","Arabidopsis_thaliana") & 
   dest_DB == "UCSC" &
   ! "UCSC" %in% names(genomeStyles("Arabidopsis thaliana"))) dest_DB <- "TAIR9"

cat("## reading inputs\n")
a <- Sys.time()
print(a)
input_DT <- as.data.table(fread(paste("zgrep","-v","^#",input_file),
                                sep="\t",
                                colClasses=c("character"),
                                header=FALSE))
print(Sys.time() - a)

print(paste("## getting names for",dest_DB))
new_names <- mapSeqlevels(unique(input_DT$V1),dest_DB)
print(Sys.time() - a)

#print(new_names)
if (length(new_names) == 0){
  print(paste("## copying data into",output_file))
  system2("cp",c(input_file,output_file))
  print(Sys.time() - a)
} else {
  print(paste("## replacing new names in DT"))
  input_DT$V1 <- new_names[input_DT$V1]
  print(Sys.time() - a)
  print(paste("## saving data with replaced names into",output_file))
  write.table(input_DT[!is.na(V1)], output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(input_DT[is.na(V1)], paste0(output_file,".chr_NA"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  print(Sys.time() - a)
}
