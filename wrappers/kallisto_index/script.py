######################################
# wrapper for rule: kallisto_index_creation
######################################
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: kallisto_index_creation \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "kallisto index -i " + snakemake.output.gen + " " + snakemake.input.cds + " >> " + snakemake.log.run  
f = open(snakemake.log.run, 'wt')
f.write("## COMMAND:\n"+version+"\n")
f.close()
shell(command)
