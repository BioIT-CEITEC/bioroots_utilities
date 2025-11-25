###########################################
# wrapper for rule: postqc_RNA_preparation
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

print("\n##\n## RULE: postqc_RNA_preparation \n##\n")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

command = "mkdir -p " + dirname(snakemake.params.bowtie2_indexes_fasta)
f.write("## COMMAND: " + command + "\n")
shell(command)

command = GTF_TO_GENPRED+" -allErrors -genePredExt "+snakemake.input.ref+" "+snakemake.output.tmp_flat+" >> "+snakemake.log.run+" 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

command = "paste <( cut -f 12 "+snakemake.output.tmp_flat+") <( cut -f 1-10 "+snakemake.output.tmp_flat+") > "+snakemake.output.flat+" 2>> "+snakemake.log.run
f.write("## COMMAND: "+command+"\n")
shell(command)

# Convert Gtf to genePred
command = GTF_TO_GENPRED+" -allErrors "+snakemake.input.ref+" "+snakemake.input.ref.replace('.gtf', '')+".genePred >> "+snakemake.log.run+" 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

# Convert genePred to bed12
command = GENE_PRED_TO_BED+" "+snakemake.input.ref.replace('.gtf','')+".genePred "+snakemake.input.ref.replace('.gtf','')+".bed12 >> "+snakemake.log.run+" 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

# sort bed12
command = "sort -k1,1 -k2,2n "+snakemake.input.ref.replace(".gtf","")+".bed12 > "+snakemake.output.bed12+" 2>> "+snakemake.log.run+" "
f.write("## COMMAND: "+command+"\n")
shell(command)

# extract rRNA data and build BOWTIE2 index
command = "cat "+ snakemake.input.ref + " | grep 'gene_biotype \"rRNA\"'  > " + snakemake.params.rRNA_prefix + ".gtf 2>> " + snakemake.log.run + " || echo '## INFO: Command returned non-zero status. Probably, there are no gene_biotype \"rRNA\" lines.' >> " + snakemake.log.run + " 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

if sum(1 for line in open(snakemake.params.rRNA_prefix + ".gtf")) == 0:
  no_rrna = True
  f.write("## INFO: file "+snakemake.params.rRNA_prefix + ".gtf is empty, therefore, skipping building of BOWTIE2 index for rRNAs."+"\n")
else:
  no_rrna = False

  #command = "sed 's/>ref|\\([^|]\\+\\)|/>\\1/' " +snakemake.input.fasta + " >> " + snakemake.log.run + " 2>&1"
  #f.write("## COMMAND: "+command+"\n")
  #shell(command)

  command = GFF_READ + " " + snakemake.params.rRNA_prefix + ".gtf -g " +snakemake.input.fasta + " -w " + snakemake.params.rRNA_prefix + ".fasta 2>> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)

  command = BOWTIE2_BUILD + " --threads " + str(snakemake.threads) + " " + snakemake.params.rRNA_prefix + ".fasta " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.params.rRNA_prefix) + ".fasta >> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)


# extract tRNA data and build BOWTIE2 index
command = "cat "+ snakemake.input.ref+" | grep 'gene_biotype \"tRNA\"' > " + snakemake.params.tRNA_prefix + ".gtf 2>> " + snakemake.log.run + " || echo '## INFO: Command returned non-zero status. Probably, there are no gene_biotype \"tRNA\" lines.' >> " + snakemake.log.run + " 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

command = "cat "+ snakemake.input.ref+" | grep 'gene_biotype \"Mt_tRNA\"' >> " + snakemake.params.tRNA_prefix + ".gtf 2>> " + snakemake.log.run + " || echo '## INFO: Command returned non-zero status. Probably, there are no gene_biotype \"Mt_tRNA\" lines.' >> " + snakemake.log.run + " 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

if sum(1 for line in open(snakemake.params.tRNA_prefix + ".gtf")) == 0:
  no_trna = True
  f.write("## INFO: file "+snakemake.params.tRNA_prefix + ".gtf is empty, therefore, skipping building of BOWTIE2 index for tRNAs."+"\n")
else:
  no_trna = False
  command = GFF_READ + " " + snakemake.params.tRNA_prefix + ".gtf -g " +snakemake.input.fasta + " -w " + snakemake.params.tRNA_prefix + ".fasta 2>> " + snakemake.log.run
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
  command = BOWTIE2_BUILD +" --threads "+ str(snakemake.threads) + " " +snakemake.input.fasta + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.input.fasta) + " >> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)
  
  # create fastq_screen.conf file
  command = "echo 'THREADS " + str(snakemake.threads) + "' > " + snakemake.output.fs_conf + " 2>> " + snakemake.log.run
  f.write("## COMMAND: "+command+"\n")
  shell(command)

  command = "echo 'DATABASE " + snakemake.params.species + " " + snakemake.params.bowtie2_indexes_fasta + os.path.basename(snakemake.input.fasta) + "' >> " + snakemake.output.fs_conf + " 2>> " + snakemake.log.run
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
