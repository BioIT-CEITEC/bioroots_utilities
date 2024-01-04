######################################
# wrapper for rule: gene_table_creation
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

command = "awk -F'[=;]' '/^>/ {{transcript=\"\"; gene=\"\"; for(i=1; i<NF; i+=2) {{ if($i==\"transcript_id\") transcript=$(i+1); if($i==\"gene_id\") gene=$(i+1); }} print transcript \"\\t\" gene }}' " + snakemake.input.cds + " > " + snakemake.output.transc
f = open(snakemake.log.run, 'wt')
f.write("\n##\n## COMMAND: " + command + "\n")
f.close()
shell(command)

command = "echo -e \"TXNAME\tGENEID\" > bob"
f = open(snakemake.log.run, 'wt')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)

command = "cat bob " + snakemake.output.transc + " > bib ; mv bib " + snakemake.output.transc + " ; rm bob"
f = open(snakemake.log.run, 'wt')
f.write("## COMMAND:\n"+command+"\n")
f.close()
shell(command)
