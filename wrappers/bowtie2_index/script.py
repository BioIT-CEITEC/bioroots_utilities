######################################
# wrapper for rule: bowtie2_index
#####################################
import os
import sys
import math
import subprocess
from os.path import basename
from os.path import dirname
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: bowtie2_index \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "bowtie2-build --threads " + str(snakemake.threads) + " -f " +snakemake.input.fasta + " " + dirname(snakemake.output.index) + "/" + basename(snakemake.input.fasta).replace(".fa", "") +  " >> " + str(snakemake.log.run) + " 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)