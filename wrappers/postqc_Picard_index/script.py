###########################################
# wrapper for rule: postqc_Picard_index
###########################################
import os
import sys
import math
import subprocess
from os.path import dirname
from snakemake.shell import shell

GTF_TO_GENPRED="gtfToGenePred"
GENE_PRED_TO_BED="genePredToBed"
GFF_READ="gffread"
# BOWTIE2_BUILD="bowtie2-build"
LOG_RUN=str(snakemake.log.run)

shell.executable("/bin/bash")

sys.stdout = open(LOG_RUN, 'a+')
f = sys.stdout

print("\n##\n## RULE: postqc_Picard_index \n##\n")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA:\n"+version+"\n")

command = GTF_TO_GENPRED+" -allErrors -genePredExt "+snakemake.input.ref+" "+snakemake.output.tmp_flat+" >> "+ LOG_RUN +" 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

command = "paste <( cut -f 12 "+snakemake.output.tmp_flat+") <( cut -f 1-10 "+snakemake.output.tmp_flat+") > "+snakemake.output.flat+" 2>> "+ LOG_RUN
f.write("## COMMAND: "+command+"\n")
shell(command)

# Convert Gtf to genePred
command = GTF_TO_GENPRED+" -allErrors "+snakemake.input.ref+" "+snakemake.input.ref.replace('.gtf', '')+".genePred >> "+ LOG_RUN +" 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

# Convert genePred to bed12
command = GENE_PRED_TO_BED+" "+snakemake.input.ref.replace('.gtf','')+".genePred "+snakemake.input.ref.replace('.gtf','')+".bed12 >> "+ LOG_RUN +" 2>&1"
f.write("## COMMAND: "+command+"\n")
shell(command)

# sort bed12
command = "sort -k1,1 -k2,2n "+snakemake.input.ref.replace(".gtf","")+".bed12 > "+snakemake.output.bed12+" 2>> "+ LOG_RUN +" "
f.write("## COMMAND: "+command+"\n")
shell(command)

