######################################
# wrapper for rule: STAR_gen_index
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

print("\n##\n## RULE: STAR_rrna_gen_index \n##\n")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

help = str(subprocess.Popen("grep -v '>' " + str(snakemake.input.gen) + " | wc -m",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
STAR_GENOME_BASES_LOG = min(14,math.floor(math.log(float(int(help)),2)/2-1))

command = "mkdir -p "+ snakemake.params.dir + " >> " + str(snakemake.log.run) + " 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)

command = "grep -Pe '\tgene\t' "+ str(snakemake.input.ref) + """ | grep -e 'rRNA' -e 'snoRNA' -e 'snRNA' | cut -d';' -f1 | sed -e 's/gene_id //' | awk 'BEGIN{FS=OFS="\t"}{print $1,$4,$5,$9,"1",$7}' > """+ str(snakemake.params.rrna_bed) + " >> " + str(snakemake.log.run) + " 2>&1
print("## COMMAND: "+command+"\n")
shell(command)

command = "bedtools getfasta -name -s -fi "+ str(snakemake.input.gen) + " -bed "+ str(snakemake.params.rrna_bed) + " -fo "+ str(snakemake.params.rrna_fa) + " >> " + str(snakemake.log.run) + " 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)

command = "STAR --runMode genomeGenerate --runThreadN "+str(snakemake.threads)+" --genomeDir "+ snakemake.params.dir +" --genomeFastaFiles "+str(snakemake.params.rrna_fa)+" --genomeSAindexNbases "+str(STAR_GENOME_BASES_LOG)+" >> "+str(snakemake.log.run)+" 2>&1 "
print("## COMMAND: "+command+"\n")
shell(command)

command = "rm -rf _STARtmp " + " >> " + str(snakemake.log.run) + " 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)

command = "cat "+ snakemake.params.log +" >> "+ str(snakemake.log.run) +" 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)
