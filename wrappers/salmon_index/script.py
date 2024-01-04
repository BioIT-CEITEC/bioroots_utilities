######################################
# wrapper for rule: salmon_index_creation
######################################
import os
import sys
import math
import subprocess
from os.path import dirname
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: salmon_index_creation \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "$(which time) --verbose sh " + snakemake.params.script + " -j 2 -a " + snakemake.input.gtf + \
          " -g " + snakemake.input.gen + " -t " + snakemake.input.cds + " -o " + dirname(snakemake.output.gen) + " >> " + snakemake.log.run
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()
shell(command)

command = "$(which time) --verbose salmon index -t " + snakemake.output.gen + " -d " + snakemake.output.dec + " -i " + snakemake.params.folder + " >> " + snakemake.log.run
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()
shell(command)
