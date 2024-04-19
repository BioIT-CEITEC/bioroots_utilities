######################################
# wrapper for rule: create_genome_dict
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: create_genome_dict \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "$(which time) picard CreateSequenceDictionary R="+snakemake.input.fasta+" O="+snakemake.output.dict+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)