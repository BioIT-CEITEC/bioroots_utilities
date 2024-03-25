######################################
# wrapper for rule: BWA_gen_index
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

sys.stdout = open(snakemake.log.run, 'a+')

print("\n##\n## RULE: convert_gff_to_gtf \n##\n")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

command = "mkdir -p "+snakemake.params.dir+" >> "+snakemake.log.run+" 2>&1"
print("## COMMAND: "+command+"\n")
shell(command)

prefix = os.path.join(snakemake.params.dir, snakemake.wildcards.ref)
command = "$(which time) bwa index -a bwtsw "+snakemake.params.extra+" -p "+prefix+" "+snakemake.input.gen+" >> "+snakemake.log.run+" 2>&1 "
print("## COMMAND: "+command+"\n")
shell(command)