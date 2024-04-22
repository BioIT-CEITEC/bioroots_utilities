######################################
# wrapper for rule: RSEM_prep_ref
######################################
import os
import sys
import math
import subprocess
# sys.path.append(os.path.abspath(os.path.dirname(__file__)+"/../../"))
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: RSEM_prep_ref \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "$(which time) rsem-prepare-reference"+\
          " --num-threads "+str(snakemake.threads)+\
          " "+snakemake.params.use_ref+" "+snakemake.input.ref+\
          " "+snakemake.input.gen+\
          " "+snakemake.params.prefix+\
          " "+snakemake.params.extra+\
          " >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)