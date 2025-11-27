###########################################
# wrapper for rule: postqc_fastq_screen_index
###########################################
import os
import sys
import math
import subprocess
from os.path import dirname
from snakemake.shell import shell

GTF_TO_GENPRED="gtfToGenePred"
GENE_PRED_TO_BED="genePredToBed"
GFF_READ="gffread"
BOWTIE2_BUILD="bowtie2-build"
LOG_RUN=str(snakemake.log.run)

# rRNA_prefix = str(snakemake.params.rRNA_prefix[0])
# tRNA_prefix = str(snakemake.params.tRNA_prefix[0])

shell.executable("/bin/bash")

# sys.stdout = open(LOG_RUN, 'a+')
# f = sys.stdout

print("\n##\n## RULE: postqc_fastq_screen_index \n##\n")
f = open(LOG_RUN, 'wt')
f.write("## RULE: postqc_fastq_screen_index \n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

command = "mkdir -p " + dirname(snakemake.params.bowtie2_indexes_fasta)
f = open(LOG_RUN, 'a+')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)

# extract rRNA data and build BOWTIE2 index
# command = "cat "+ snakemake.input.ncbi_annot + " | grep 'gbkey=rRNA' | grep -v 'ribosomal RNA protein' > " + snakemake.params.rRNA_prefix.replace("fasta","gff") + " 2>> " + LOG_RUN + " || echo '## INFO: Command returned non-zero status. Probably, there are no gbkey=rRNA lines.' >> " + LOG_RUN + " 2>&1"
command = "cat "+ snakemake.input.ncbi_annot + " | grep 'gbkey=rRNA' | grep -v 'ribosomal RNA protein' > " + snakemake.params.rRNA_prefix.replace("fasta","gff")
f = open(LOG_RUN, 'a+')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)

if sum(1 for line in open(snakemake.params.rRNA_prefix.replace("fasta","gff"))) == 0:
  no_rrna = True
  f = open(LOG_RUN, 'a+')
  f.write("## INFO: file "+ snakemake.params.rRNA_prefix.replace("fasta","gff") + " is empty, therefore, skipping building of BOWTIE2 index for rRNAs."+"\n")
  f.close()
else:
  no_rrna = False
  f = open(LOG_RUN, 'a+')
  f.write("## INFO: file "+ snakemake.params.rRNA_prefix.replace("fasta","gff") + " is not empty, proceeding with building of BOWTIE2 index for rRNAs."+"\n")
  f.close()

  #command = "sed 's/>ref|\\([^|]\\+\\)|/>\\1/' " +snakemake.input.ncbi_genomic + " >> " + LOG_RUN + " 2>&1"
  #f.write("## COMMAND: "+command+"\n")
  #shell(command)

  command = GFF_READ + " " + str(snakemake.params.rRNA_prefix).replace("fasta","gff") + " -g " + str(snakemake.input.ncbi_genomic) + " -w " + str(snakemake.params.rRNA_prefix) + " 2>> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)

  command = BOWTIE2_BUILD + " --threads " + str(snakemake.threads) + " " + rRNA_prefix + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.rRNA_prefix) + " >> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)


# extract tRNA data and build BOWTIE2 index
command = "cat "+ snakemake.input.ncbi_annot+" | grep 'gbkey=tRNA' > " + snakemake.params.tRNA_prefix.replace("fasta","gff") + " 2>> " + LOG_RUN + " || echo '## INFO: Command returned non-zero status. Probably, there are no gbkey=tRNA lines.' >> " + LOG_RUN + " 2>&1"
f = open(LOG_RUN, 'a+')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)

if sum(1 for line in open(snakemake.params.tRNA_prefix.replace("fasta","gff"))) == 0:
  no_trna = True
  f.write("## INFO: file "+snakemake.params.tRNA_prefix.replace("fasta","gff") + " is empty, therefore, skipping building of BOWTIE2 index for tRNAs."+"\n")
else:
  no_trna = False
  command = GFF_READ + " " + str(snakemake.params.tRNA_prefix).replace("fasta","gff") + " -g " + str(snakemake.input.ncbi_genomic) + " -w " + str(snakemake.params.tRNA_prefix) + " 2>> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)

  command = BOWTIE2_BUILD +" --threads "+ str(snakemake.threads) + " " + snakemake.params.tRNA_prefix + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.tRNA_prefix) + " >> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)


if no_rrna and no_trna:
  # there are no tRNA nor rRNA sequences so an empty fastq_screen.conf is generated
  f = open(LOG_RUN, 'a+')
  f.write("## INFO: there are no tRNA nor rRNA sequences so an empty fastq_screen.conf is generated\n")
  f.close()
    
  command = "touch "+ snakemake.output.fs_conf + " >> " + LOG_RUN + " 2>&1"
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)
else:
  # build BOWTIE2 index for whole genome
  command = BOWTIE2_BUILD +" --threads "+ str(snakemake.threads) + " " + snakemake.input.ncbi_genomic + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.input.ncbi_genomic) + " >> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)
  
  # create fastq_screen.conf file
  command = "echo 'THREADS " + str(snakemake.threads) + "' > " + snakemake.output.fs_conf + " 2>> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)

  command = "echo 'DATABASE " + snakemake.params.species + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.input.ncbi_genomic) + "' >> " + snakemake.output.fs_conf + " 2>> " + LOG_RUN
  f = open(LOG_RUN, 'a+')
  f.write("## COMMAND:\n"+command+"\n")
  f.close()
  shell(command)

  if not no_rrna:
    command = "echo 'DATABASE rRNA " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.rRNA_prefix) + "' >> " + snakemake.output.fs_conf + " 2>> " + LOG_RUN
    f = open(LOG_RUN, 'a+')
    f.write("## COMMAND:\n"+command+"\n")
    f.close()
    shell(command)
    
  if not no_trna:
    command = "echo 'DATABASE tRNA " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.tRNA_prefix) + "' >> " + snakemake.output.fs_conf + " 2>> " + LOG_RUN
    f = open(LOG_RUN, 'a+')
    f.write("## COMMAND:\n"+command+"\n")
    f.close()
    shell(command)
