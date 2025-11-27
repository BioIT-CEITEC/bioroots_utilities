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

shell.executable("/bin/bash")

sys.stdout = open(snakemake.log.run, 'a+')
f = sys.stdout

print("\n##\n## RULE: postqc_fastq_screen_index \n##\n")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

command = "mkdir -p " + dirname(snakemake.params.bowtie2_indexes_fasta)
f.write("## COMMAND: " + command + "\n")
shell(command)

# extract rRNA data and build BOWTIE2 index
# command = "cat "+ snakemake.input.ncbi_annot + " | grep 'gbkey=rRNA' | grep -v 'ribosomal RNA protein' > " + snakemake.params.rRNA_prefix + ".gff 2>> " + snakemake.log.run + " || echo '## INFO: Command returned non-zero status. Probably, there are no gbkey=rRNA lines.' >> " + snakemake.log.run + " 2>&1"
command = "cat "+ snakemake.input.ncbi_annot + " | grep 'gbkey=rRNA' | grep -v 'ribosomal RNA protein' > " + snakemake.params.rRNA_prefix + ".gff
f.write("## COMMAND: "+command+"\n")
shell(command)

if sum(1 for line in open(snakemake.params.rRNA_prefix + ".gff")) == 0:
  no_rrna = True
  f.write("## INFO: file "+snakemake.params.rRNA_prefix + ".gff is empty, therefore, skipping building of BOWTIE2 index for rRNAs."+"\n")
else:
  no_rrna = False

  #command = "sed 's/>ref|\\([^|]\\+\\)|/>\\1/' " +snakemake.input.ncbi_genomic + " >> " + snakemake.log.run + " 2>&1"
  #f.write("## COMMAND: "+command+"\n")
  #shell(command)

  command = GFF_READ + " " + snakemake.params.rRNA_prefix + ".gff -g " +snakemake.input.ncbi_genomic + " -w " + snakemake.params.rRNA_prefix + ".fasta 2>> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)

  command = BOWTIE2_BUILD + " --threads " + str(snakemake.threads) + " " + snakemake.params.rRNA_prefix + ".fasta " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.rRNA_prefix) + ".fasta >> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)


# extract tRNA data and build BOWTIE2 index
command = "cat "+ snakemake.input.ncbi_annot+" | grep 'gbkey=tRNA' > " + snakemake.params.tRNA_prefix + ".gff 2>> " + snakemake.log.run + " || echo '## INFO: Command returned non-zero status. Probably, there are no gbkey=tRNA lines.' >> " + snakemake.log.run + " 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

if sum(1 for line in open(snakemake.params.tRNA_prefix + ".gff")) == 0:
  no_trna = True
  f.write("## INFO: file "+snakemake.params.tRNA_prefix + ".gff is empty, therefore, skipping building of BOWTIE2 index for tRNAs."+"\n")
else:
  no_trna = False
  command = GFF_READ + " " + snakemake.params.tRNA_prefix + ".gff -g " +snakemake.input.ncbi_genomic + " -w " + snakemake.params.tRNA_prefix + ".fasta 2>> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)

  command = BOWTIE2_BUILD +" --threads "+ str(snakemake.threads) + " " + snakemake.params.tRNA_prefix + ".fasta " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.tRNA_prefix) + ".fasta >> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)


if no_rrna and no_trna:
  # there are no tRNA nor rRNA sequences so an empty fastq_screen.conf is generated
  f.write("## INFO: there are no tRNA nor rRNA sequences so an empty fastq_screen.conf is generated\n")
  
  command = "touch "+ snakemake.output.fs_conf + " >> " + snakemake.log.run + " 2>&1"
  f.write("## COMMAND: "+command+"\n")
  shell(command)
else:
  # build BOWTIE2 index for whole genome
  command = BOWTIE2_BUILD +" --threads "+ str(snakemake.threads) + " " +snakemake.input.ncbi_genomic + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.input.ncbi_genomic) + " >> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)
  
  # create fastq_screen.conf file
  command = "echo 'THREADS " + str(snakemake.threads) + "' > " + snakemake.output.fs_conf + " 2>> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)

  command = "echo 'DATABASE " + snakemake.params.species + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.input.ncbi_genomic) + "' >> " + snakemake.output.fs_conf + " 2>> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)

  if not no_rrna:
    command = "echo 'DATABASE rRNA " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.rRNA_prefix) + ".fasta' >> " + snakemake.output.fs_conf + " 2>> " + snakemake.log.run
    f.write("## COMMAND: "+command+"\n")
    shell(command)
    
  if not no_trna:
    command = "echo 'DATABASE tRNA " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.tRNA_prefix) + ".fasta' >> " + snakemake.output.fs_conf + " 2>> " + snakemake.log.run
    f.write("## COMMAND: "+command+"\n")
    shell(command)
