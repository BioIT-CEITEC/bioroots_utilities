######################################
# wrapper for rule: gtf_to_fasta
######################################
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'wt')
f.write("\n##\n## RULE: gtf_to_fasta \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "gffread -F " + snakemake.input.gtf + " -g " + snakemake.input.gen + " -x " + snakemake.output.cds
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()
shell(command)

command = "gffread -F " + snakemake.input.gtf + " -g " + snakemake.input.gen + " -w " + snakemake.output.cdna
f = open(snakemake.log.run, 'wt')
f.write("## CONDA:\n"+version+"\n")
f.close()
shell(command)
