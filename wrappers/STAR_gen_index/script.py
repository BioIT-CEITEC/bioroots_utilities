######################################
# wrapper for rule: STAR_gen_index
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

print("\n##\n## RULE: STAR_gen_index \n##\n")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

help = str(subprocess.Popen("grep -v '>' " + str(snakemake.input.gen) + " | wc -m",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
STAR_GENOME_BASES_LOG = min(14,math.floor(math.log(float(int(help)),2)/2-1))

command = "mkdir -p "+ snakemake.params.dir + " >> " + str(snakemake.log.run) + " 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)

command = "STAR --runMode genomeGenerate --runThreadN "+str(snakemake.threads)+" --limitGenomeGenerateRAM " + str(snakemake.resources.mem * 1000000000) +" --genomeDir "+ snakemake.params.dir +" --genomeFastaFiles "+str(snakemake.input.gen)+" --sjdbGTFfile " + str(snakemake.input.ref) + " --genomeSAindexNbases "+str(STAR_GENOME_BASES_LOG)+" "+str(snakemake.params.extra)+" >> "+str(snakemake.log.run)+" 2>&1 "
print("## COMMAND: "+command+"\n")
shell(command)

command = "rm -rf _STARtmp " + " >> " + str(snakemake.log.run) + " 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)

command = "cat "+ snakemake.params.log +" >> "+ str(snakemake.log.run) +" 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)
