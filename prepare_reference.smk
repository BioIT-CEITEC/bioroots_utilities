import math
import subprocess
import json
import re
import os.path
from snakemake.utils import R
from snakemake.utils import report
from os.path import split



wildcard_constraints:
    donor="[^\.]+",\
    sampling="[^\.]+",\
    analysis="[^\.]+",\
    tag="[^\.]+",\
    organism="[^\.]+",\
    reference="[^\.]+",\
    germline="[^\.]+",\
    somatic="[^\.]+",\
    owner="[^\.]+",\
    run_name="[^\.]+",\
    sample="[^\.]+",\
    tech_rep="[^\.]+",\
    tech_rep_isnone="[^\.]+"


########################################
# DEFINITION OF usefull FUNCTIONS
#


########################################
# DEFINITION OF CONFIGURATION FEATURES
#


####################################
# DEFINITION OF TOOLS
#
DELLY = "delly"
FASTQC = "fastqc"
GLOBAL_ADAPTERS = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/adapters_merge.fa"
REAPER_SRC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/reaper-15-065/src"
TRIMMOMATIC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/trimmomatic-master/classes/trimmomatic.jar"
GFFREAD = "gffread"
STAR = "STAR"
SAMTOOLS = "samtools"
RSEM_PATH = ""
BBMAP = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/bbmap"
FEATURE_COUNTS = "featureCounts"
UCSC_SCRIPTS = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools"
GTF_TO_BED12 = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/gtf2bed12.py"
PICARD = "picard"
PRESEQ = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/preseq_v2.0.2/preseq" #"/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/preseq_v2.0.1/preseq"
RSEQC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/RSeQC-2.6.4/build/scripts-2.7"
DUPRADAR = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/dupRadar.R"
# New and needed
convert_to_ucsc = workflow.basedir + "/../scripts/convert_chromosome_names.R"


####################################
# DEFINITION OF FINAL RULES
#

# rule interval_list_prep:
#     input:  bed = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.bed",
#             bam =
#     output: idx = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.interval_list"
#     run:

#rule for picard bed index selection
#for DNA





# OTHER PREP RULES
#
#rule snp_vcf2bed:
#    input:  vcf = "{dir}/{species}/{ref}/other/snp/{ref}.snp.vcf"
#    output: bed = "{dir}/{species}/{ref}/other/snp/{ref}.snp.bed"
#    params: sort_dir = "/mnt/ssd/ssd_1/tmp/"
#    resources:  mem = 10
#    shell:
#        "convert2bed -i vcf -m {resources.mem}G -r {params.sort_dir} < {input.vcf} > {output.bed}"


#rule makeblastdb_for_crispr_csv:
#    input:  csv = "{dir}/{species}/{ref}/{ref}.csv",
#    output: fa  = "{dir}/{species}/{ref}/index/{ref}.uniq_revcomp_ins.fa",
#            mod = "{dir}/{species}/{ref}/{ref}_mod.csv",
#            uniq= "{dir}/{species}/{ref}/{ref}.uniq_revcomp_ins.csv",
#    params: script = workflow.basedir+"/../wraps/prepare_reference/makeblastdb_for_crispr_csv/uniq_revcomp_crispr.R",
#    log:    run = "{dir}/{species}/{ref}/{ref}.makeblastdb_for_crispr_csv.log",
#    conda:  "../wraps/prepare_reference/makeblastdb_for_crispr_csv/env.yaml"
#    script: "../wraps/prepare_reference/makeblastdb_for_crispr_csv/script.py"


# rule biotype_groups_download:
#     input:  info = "{dir}/{species}/{ref}/info.txt"
#     output: bio_groups = "{dir}/{species}/{ref}/annot/{ref}.biotype_groups.txt"
#     run:
#         file = open(input.info, "r")
#         for line in file:
#              if re.search("release::", line):
#                  version = line.replace("release::","")
#
#         shell("mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_{version} -e \"select distinct(name),biotype_group from biotype where db_type like '%core%' order by biotype_group,name;\" > {output.bio_groups}")
#

#rule gtf_to_ensembl_DB:
#    input:  gtf = "{dir}/{organism}/{ref}/annot/{ref}.gtf",
#    output: sqlite = "{dir}/{organism}/{ref}/annot/{ref}.sqlite.gz",
#    params: dir = workflow.basedir,
#            info = "{dir}/{organism}/{ref}/info.txt"
#    conda:  "../wraps/prepare_reference/gtf_to_ensembl_DB/env.yaml"
#    shell:
#        "Rscript {params.dir}/../wraps/prepare_reference/gtf_to_ensembl_DB/script.R {wildcards.organism} {wildcards.ref} {params.info} {input.gtf} {output.sqlite}"

# KIT INTERVALS RULES
#
#rule kit_snp_list:
#    input:  bed = "{dir}/{organism}/{ref}/intervals/{kit}/{kit}.bed",
#            snp_bed = "{dir}/{organism}/{ref}/other/snp/{ref}_1k_gen_0.1_snps_pos.bed",
#    output: kit_snp_bed = "{dir}/{organism}/{ref}/other/snp/{kit}/{kit}_snps.bed",
#    log:    run = "{dir}/{organism}/{ref}/other/snp/{kit}/kit_snp_list.log",
#    conda:  "../wraps/prepare_reference/kit_snp_list/env.yaml"
#    script: "../wraps/prepare_reference/kit_snp_list/script.py"


#rule index_and_gz_bed:
#    input:  bed = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.bed",
#    output: bed_gz = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.bed.gz",
#            bed_gz_tbi = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.bed.gz.tbi",
#    run:
#        shell("bgzip -c {input.bed} > {output.bed_gz}")
#        shell("tabix -p bed {output.bed_gz}")

#rule interval_list_from_bed:
#    input:  bed = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.bed",
#            dict = "{dir}/{species}/{ref}/seq/{ref}.dict",
#    output: interval_list = "{dir}/{species}/{ref}/intervals/{kit}/{kit}.interval_list",
#    log:    run = "{dir}/{species}/{ref}/intervals/{kit}/Picard_BedToIntervalList_run.log",
#    conda:  "../wraps/prepare_reference/interval_list_from_bed/env.yaml"
#    script: "../wraps/prepare_reference/interval_list_from_bed/script.py"

#rule bed_from_dict_for_wgs:
#    input:  dict = "{dir}/{species}/{ref}/seq/{ref}.dict",
#    output: bed = "{dir}/{species}/{ref}/intervals/wgs/wgs.bed",
#    conda:  "../wraps/prepare_reference/bed_from_dict_for_wgs/env.yaml"
#    script: "../wraps/prepare_reference/bed_from_dict_for_wgs/script.py"

#rule create_genome_dict:
#    input:  fasta = "{dir}/{species}/{ref}/seq/{ref}.fa",
#    output: dict = "{dir}/{species}/{ref}/seq/{ref}.dict",
#    log:    run = "{dir}/{species}/{ref}/seq/create_genome_dict_run.log",
#    conda:  "../wraps/prepare_reference/create_genome_dict/env.yaml"
#    script: "../wraps/prepare_reference/create_genome_dict/script.py"


#rule create_IGV_genome:
#    input:  fa = "{dir}/{species}/{ref}/seq/{ref}.fa",
#            info = "{dir}/{species}/{ref}/info.txt",
#    output: prop = "{dir}/{species}/{ref}/seq/property.txt",
#            genome = "{dir}/{species}/{ref}/seq/{ref}.genome",
#            # html = "{dir}/{species}/{ref}/seq/{ref}.gen_for_IGV.html",
#    log:    run = "{dir}/{species}/{ref}/seq/{ref}.gen_for_IGV.log",
#    params: gtf = "{dir}/{species}/{ref}/annot/{ref}.gtf",
#            gff = "{dir}/{species}/{ref}/annot/{ref}.gff3",
#            dir = "{dir}/{species}/{ref}/seq/",
#    conda:  "../wraps/prepare_reference/create_IGV_genome/env.yaml"
#    script: "../wraps/prepare_reference/create_IGV_genome/script.py"


#rule download_VEP_data:
#    input:  info = "{dir}/{species}/{ref}/info.txt",
#    output: html = "{dir}/{species}/{ref}/annot/vep/{ref}.download_vep_data.html",
#    log:    run =  "{dir}/{species}/{ref}/annot/vep/{ref}.download_vep_data.log",
#    params: dir =  "{dir}/{species}/{ref}/annot/vep/",
#            ver =  ""
#    conda:  "../wraps/prepare_reference/download_VEP_data/env.yaml"
#    script: "../wraps/prepare_reference/download_VEP_data/script.py"
#rule download_VEP_data_merged:
#    input:  info = "{dir}/{species}/{ref}/info.txt",
#    output: html = "{dir}/{species}/{ref}/annot/vep/{ref}.download_vep_data_merged.html",
#    log:    run =  "{dir}/{species}/{ref}/annot/vep/{ref}.download_vep_data_merged.log",
#    params: dir =  "{dir}/{species}/{ref}/annot/vep/",
#            ver =  "merged"
#    conda:  "../wraps/prepare_reference/download_VEP_data/env.yaml"
#    script: "../wraps/prepare_reference/download_VEP_data/script.py"
#rule download_VEP_data_refseq:
#    input:  info = "{dir}/{species}/{ref}/info.txt",
#    output: html = "{dir}/{species}/{ref}/annot/vep/{ref}.download_vep_data_refseq.html",
#    log:    run =  "{dir}/{species}/{ref}/annot/vep/{ref}.download_vep_data_refseq.log",
#    params: dir =  "{dir}/{species}/{ref}/annot/vep/",
#            ver =  "refseq"
#    conda:  "../wraps/prepare_reference/download_VEP_data/env.yaml"
#    script: "../wraps/prepare_reference/download_VEP_data/script.py"

#rule process_gtf_for_VEP:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gtf"
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.vep_ready.gtf.gz",
#            tbi = "{dir}/{species}/{ref}/annot/{ref}.vep_ready.gtf.gz.tbi"
#    shell:
#        """
#        grep -v "#" {input.ref} | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output.ref}
#        tabix -p gff {output.ref}
#        """

#rule process_gff_for_VEP:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gff3"
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.vep_ready.gff3.gz",
#            tbi = "{dir}/{species}/{ref}/annot/{ref}.vep_ready.gff3.gz.tbi"
#    shell:
#        """
#        grep -v "#" {input.ref} | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output.ref}
#        tabix -p gff {output.ref}
#        """

#rule process_transdecoder_ids_in_gtf:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gtf"
#    output: tab = "{dir}/{species}/{ref}/annot/{ref}.transdecoder_ids_map"
#    log:    run = "{dir}/{species}/{ref}/annot/{ref}.transdecoder_ids_map.log"
#    params: rscript = workflow.basedir+"/../wraps/prepare_reference/process_transdecoder_ids_in_gtf/process_gtf_ids.R"
#    conda:  "../wraps/prepare_reference/process_transdecoder_ids_in_gtf/env.yaml"
#    script: "../wraps/prepare_reference/process_transdecoder_ids_in_gtf/script.py"



# GENOME INDEXING RULE
# (if genome index is not defined in config.json it will be generated from the scratch,
# otherwise, it will be just soft-linked from the source)

#rule smallRNA_prep_contam:
#    input:  ref = "{dir}/{species}/smallRNA/contaminants/{rna_type}_contaminants.fasta",
#    output: counts   = "{dir}/{species}/smallRNA/contaminants/{rna_type}_contaminants.fasta.counts",
#            bin_counts = "{dir}/{species}/smallRNA/contaminants/{rna_type}_contaminants.fasta.counts.obinary",
#            header   = "{dir}/{species}/smallRNA/contaminants/{rna_type}_contaminants.fasta.header",
#            #html = "{dir}/{species}/smallRNA/contaminants/{rna_type}_contaminants.indexation_run.html",
#    log:    run = "{dir}/{species}/smallRNA/contaminants/{rna_type}_contaminants.indexation_run.log",
#    threads:    20,
#    params: extra = "",
#            dir = "{dir}/{species}/smallRNA/contaminants/",
#    conda:  "../wraps/prepare_reference/smallRNA_prep_contam/env.yaml"
#    script: "../wraps/prepare_reference/smallRNA_prep_contam/script.py"


rule postqc_RNA_preparation:
    input:  ref = config["organism_gtf"],
            ncbi_annot = config["organism_ncbi_gff"],
            ncbi_genomic = config["organism_ncbi_general"]
    output: bed12 = config["organism_picard_bed12"],
            tmp_flat = temp(config["organism_picard_refFlat"]+".tmp"),
            flat = config["organism_picard_refFlat"],
            fs_conf = config["reference_dir"]+"/seq/BOWTIE2_fastq_screen/fastq_screen.conf",
    log:    run = config["reference_dir"] + "/annot/" + config["release"] + "/Picard/Picard_preparation.log",
    threads:   15,
    params: species = config["organism"],
            bowtie2_indexes_fasta = config["reference_dir"] + "/seq/BOWTIE2_fastq_screen/",
            rRNA_prefix = config["organism_ncbi_rRNA"],
            tRNA_prefix = config["organism_ncbi_tRNA"],
    conda:  "../wrappers/postqc_RNA_preparation/env.yaml"
    script: "../wrappers/postqc_RNA_preparation/script.py"

rule BWA_gen_index:
    input:  gen = config["organism_fasta"],
            idx = config["organism_fasta"]+".fai",
    output: bwt = config["organism_bwa"],
    log:    run = config["organism_kallisto"] + "/tool_dir/BWA/BWA.indexation_run.log",
    threads:    20
    params: extra = "",
            dir = config["reference_dir"]+"/tool_dir/BWA/"
    conda:  "../wrappers/BWA_gen_index/env.yaml"
    script: "../wrappers/BWA_gen_index/script.py"

rule STAR_gen_index:
    input:  gen = config["organism_fasta"],
            idx = config["organism_fasta"]+".fai",
            ref = config["organism_gtf"],
    output: SAindex = config["organism_star"],
    params: dir = config["reference_dir"]+"/tool_data/STAR/"+config["release"],
            log = config["reference_dir"]+"/tool_data/STAR/"+config["release"]+"/Log.out",
            extra = "",
    resources:  mem = 100
    log:    run = config["reference_dir"]+"/tool_data/STAR/"+config["release"]+"/"+config["release"]+".indexation_run.log",
    threads:    30
    conda:  "../wrappers/STAR_gen_index/env.yaml"
    script: "../wrappers/STAR_gen_index/script.py"

rule create_salmon_index:
  input:  gen = config["organism_fasta"],
          gtf = config["organism_gtf"],
          cdna = config["organism_cdna_fasta"]
  output: gen = config["organism_salmon_gentrome"],
          dec = config["organism_salmon"]+"/decoy/decoys.txt",
  log:    run = config["organism_salmon"]+"/"+config["release"]+"_salmon.decoy_creation.log"
  params: folder = config["organism_salmon"]
  threads: 20
  conda: "../wrappers/salmon_index/env.yaml"
  script: "../wrappers/salmon_index/script.py"

rule create_kallisto_index:
  input:  cdna = config["organism_cdna_fasta"],
  output: gen = config["organism_kallisto"],
  log:    run = config["reference_dir"]+"/tool_data/Kallisto/" + config["release"]+ "/Kallisto.decoy_creation.log",
  threads: 20
  conda: "../wrappers/kallisto_index/env.yaml"
  script: "../wrappers/kallisto_index/script.py"

rule create_gene_table:
  input:  cdna = config["organism_cdna_fasta"],
  output: transc = config["reference_dir"] + "/annot/" + config["release"] + "/transcript_gene.txt",
  log:    run = config["reference_dir"]+"/annot/"+config["release"]+"/"+config["release"]+"_gene_table.log",
  script: "../wrappers/gene_table/script.py"

rule gtf_to_fasta:
    input:  gen = config["organism_fasta"],
            gtf = config["organism_gtf"]
    output: cds = config["organism_cds_fasta"],
            cdna = config["organism_cdna_fasta"],
    log:    run = config["reference_dir"] + "/annot/" + config["release"] + "/" + config["assembly"] + ".log"
    conda: "../wrappers/gtf_to_fasta/env.yml"
    script: "../wrappers/gtf_to_fasta/script.py"

rule bowite2_index:
    input: fasta = config["organism_fasta"]
    output: index = config["organism_bowtie2"]
    log: run = config["reference_dir"] + "/tool_data/Bowtie2/bowtie2_index_run.log"
    threads: 30
    conda: "../wrappers/bowtie2_index/env.yaml"
    script: '../wrappers/bowtie2_index/script.py'

rule RSEM_prep_ref:
    input:  ref = config["organism_gtf"],
            gen = config["organism_fasta"],
    output: idx = config["organism_rsem"],
    log:    run = config["reference_dir"] + "/tool_data/RSEM/" + config["release"] + "/" + config["release"] + ".indexation_run.log",
    threads:    10,
    params: prefix = config["reference_dir"] + "/tool_data/RSEM/" + config["release"] + "/" + config["assembly"],
            dir = config["reference_dir"] + "/tool_data/RSEM/" + config["release"],
            use_ref = "--gtf",
            extra = "",
    conda:  "../wrappers/RSEM_prep_ref/env.yaml"
    script: "../wrappers/RSEM_prep_ref/script.py"

rule make_fasta_idx_ucsc_version:
    input:  idx = config["organism_fasta"] + ".fai",
    output: ucsc = config["organism_fasta"] + ".fai.ucsc",
    log:    run = config["reference_dir"] + "/seq/" + config["assembly"] + ".make_fasta_idx_ucsc_version.log",
    params: dest = "UCSC",
            organism = config["organism"],
    conda:  "../wrappers/make_fasta_idx_ucsc_version/env.yaml"
    script: "../wrappers/make_fasta_idx_ucsc_version/script.R"

rule chrom_sizes:
    input:  idx = config["organism_fasta"] + ".fai",
    output: chs = config["organism_chr_sizes"],
    shell:  "cut -f 1,2 {input.idx} > {output.chs}"

rule create_genome_dict:
    input:  fasta = config["organism_fasta"],
    output: dict = config["organism_dict"],
    log:    run = config["reference_dir"] + "/seq/create_genome_dict_run.log",
    conda:  "../wrappers/create_genome_dict/env.yaml"
    script: "../wrappers/create_genome_dict/script.py"

rule make_fasta_idx:
    input:  gen = config["organism_fasta"],
    output: idx = config["organism_fasta"] + ".fai",
    threads:  1
    conda:  "../wrappers/STAR_gen_index/env.yaml"
    shell:
        "samtools faidx {input.gen}"

# rule make_fasta_idx:
#     input:  gen = config["organism_ncbi_general"],
#     output: idx = config["organism_ncbi_general"] + ".fai",
#     threads:  1
#     conda:  "../wrappers/STAR_gen_index/env.yaml"
#     shell:
#         "samtools faidx {input.gen}"

# COMBINE REF FOR SPIKEIN
#

#def combine_spike_in_fasta_input(wildcards):
#    if "_spike_" in wildcards.ref:
#        orig_ref_name = wildcards.ref.split("_spike_")[0]
#        orig_ref_fasta = os.path.join(wildcards.dir,wildcards.species,orig_ref_name,"seq",orig_ref_name + ".fa")
#        spike_ref_name = wildcards.ref.split("_spike_")[1]
#        spike_ref_fasta = os.path.join(wildcards.dir,"spike_ins",spike_ref_name,"seq",spike_ref_name + ".fa")
#        return { 'orig_fa': orig_ref_fasta,\
#                 'spike_fa': spike_ref_fasta}
#    else:
#        return { 'orig_fa': "not_spike_in_reference_should_exist",\
#                 'spike_fa': "not_spike_in_reference_should_exist"}

#rule combine_spike_in_fasta:
#    input:  unpack(combine_spike_in_fasta_input),
#    output: fa = "{dir}/{species}/{ref}/seq/{ref}.fa",
#    params: ref_dir = "{dir}/{species}/{ref}"
#    resources:  mem = 10
#    shell:
#        """
#        mkdir -p {params.ref_dir}/seq
#        cat {input.orig_fa} {input.spike_fa} > {output.fa}
#        """

#def combine_spike_in_gtf_input(wildcards):
#    if "_spike_" in wildcards.ref:
#        orig_ref_name = wildcards.ref.split("_spike_")[0]
#        orig_ref_gtf = os.path.join(wildcards.dir,wildcards.species,orig_ref_name,"annot",orig_ref_name + ".gtf")
#        spike_ref_name = wildcards.ref.split("_spike_")[1]
#        spike_ref_gtf = os.path.join(wildcards.dir,"spike_ins",spike_ref_name,"annot",spike_ref_name + ".gtf")
#        return { 'orig_gtf': orig_ref_gtf,\
#                 'spike_gtf': spike_ref_gtf,}
#    else:
#        return { 'orig_gtf': "not_spike_in_reference_should_exist",\
#                 'spike_gtf': "not_spike_in_reference_should_exist",}

#rule combine_spike_in_gtf:
#    input:  unpack(combine_spike_in_gtf_input),
#    output: gtf = "{dir}/{species}/{ref}/annot/{ref}.gtf",
#    params: ref_dir = "{dir}/{species}/{ref}"
#    resources:  mem = 10
#    shell:
#        """
#        mkdir -p {params.ref_dir}/annot
#        cat {input.orig_gtf} {input.spike_gtf} > {output.gtf}
#        """

# GENOMES MANAGER RULES
#
### TODO: add check between genome and annotation (e.g., chromosome names)

# rule extract_biotype_list:
#     input:  info = "{dir}/{species}/{ref}/info.txt",
#     output: list = "{dir}/{species}/{ref}/annot/biotype_list.txt",
#     run:
#
#         shell("""
#             REL=`cat {input.info} | grep 'release::' | cut -f 3 -d ':'`
#             mysql -uanonymous -P3306 -hensembldb.ensembl.org -Densembl_production_$REL -e "select distinct(name),biotype_group from biotype where db_type like '%core%' and is_current=1 order by biotype_group,name;" > {output.list}
#             # what if it's not possible to generate one
#         """)


#rule bed_ucsc_to_ensembl:
#    input:  bed = "{dir}/{species}/{ref}/other/a2i_data/ucsc/{file}.gz",
#    output: bed = "{dir}/{species}/{ref}/other/a2i_data/ensembl/{file}.gz",
#    log:    run = "{dir}/{species}/{ref}/other/a2i_data/ensembl/{file}.bed_ucsc_to_ensembl.log",
#    params: bed = "{dir}/{species}/{ref}/other/a2i_data/ensembl/{file}",
#    conda:  "../wraps/prepare_reference/make_fasta_idx_ucsc_version/env.yaml"
#    shell:
#        "Rscript "+convert_to_ucsc+" {input.bed} {params.bed} Ensembl {wildcards.species} > {log.run} 2>&1 && gzip {params.bed}"




#rule make_fasta_idx_ncbi:
#    input:  gen = "{dir}/{species}/{ref}/annot/ncbi/{ref}.ncbi.fna",
#    output: idx = "{dir}/{species}/{ref}/annot/ncbi/{ref}.ncbi.fna.fai",
#    threads:  1
#    conda:  "../wraps/prepare_reference/make_fasta_idx/env.yaml"
#    shell:
#        "samtools faidx {input.gen}"

#
# rule copy_genome:
#   input:  IN_GENOMES_LIST
#   output: OUT_GENOMES_LIST
#   threads:  1
#   run:
#     for inp, out in zip(input, output):
#         shell(" ln {inp} {out} ")


# REFERENCES MANAGER RULES
#

#rule ref_gtf_to_bed_intervals:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gtf",
#    output: ref = "{dir}/{species}/{ref}/intervals/rna/rna.bed",
#    params: rscript = workflow.basedir + "/../wraps/prepare_reference/ref_gtf_to_bed_intervals/gtf_to_bed_intervals.R"
#    threads:  1
#    conda:  "../wraps/prepare_reference/ref_gtf_to_bed_intervals/env.yaml"
#    script: "../wraps/prepare_reference/ref_gtf_to_bed_intervals/script.py"


#rule ref_gtf_to_bed:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gtf",
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.bed",
#    params: sort_dir = "/mnt/ssd/ssd_1/tmp/",
#    threads:  1
#    resources:  mem = 10  # 10 GB
#    shell:
#        "awk '$3 !~ \"gene\"' {input.ref} | convert2bed -i gtf -m {resources.mem}G -r {params.sort_dir} - > {output.ref}"

#rule ref_gff_to_bed:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gff3",
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.bed",
#    params: sort_dir = "/mnt/ssd/ssd_1/tmp/",
#    threads:  1
#    resources:  mem = 10  # 10 GB
#    shell:
#        "convert2bed -i gff -m {resources.mem}G -r {params.sort_dir} < {input.ref} > {output.ref}"

#rule unzip_ref_gtf:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gtf.gz"
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.gtf",
#    threads:  10
#    shell:
#        "unpigz -p {threads} -c {input.ref} > {output.ref} "

#rule ref_gff3_to_gtf:
#    input:  ref = ancient("{dir}/{species}/{ref}/annot/{ref}.gff3"),
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.gtf.gz",
#    params: ref = "{dir}/{species}/{ref}/annot/{ref}.temp.gtf",
#            snakemake_dir = workflow.basedir + "/../",
#    log:    run = "{dir}/{species}/{ref}/annot/{ref}_gff3_to_gtf.log",
#    conda:  "../wraps/prepare_reference/gffread_rules/env.yaml"
#    script: "../wraps/prepare_reference/gffread_rules/script.py"

# rule ref_gff_to_gtf:
#     input:  ref = ancient("{dir}/{species}/{ref}/annot/{ref}.gff"),
#     output: ref = "{dir}/{species}/{ref}/annot/{ref}.gtf",
#     params: ref = "{dir}/{species}/{ref}/annot/{ref}.temp.gtf",
#             snakemake_dir = workflow.basedir + "/../",
#     conda:  "../wraps/prepare_reference/gffread_rules/env.yaml"
#     script: "../wraps/prepare_reference/gffread_rules/script.py"
#     #shell:
#     #    "gffread {input.ref} -T -O -o - | sed -r 's/([^_])biotype/\1gene_biotype/' > {output.ref} "

#rule unzip_ref_gff3:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gff3.gz"
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.gff3",
#    threads:  10
#    shell:
#        "unpigz -p {threads} -c {input.ref} > {output.ref} "

#rule unzip_ref_gff:
#    input:  ref = "{dir}/{species}/{ref}/annot/{ref}.gff.gz"
#    output: ref = "{dir}/{species}/{ref}/annot/{ref}.gff",
#    threads:  10
#    shell:
#        "unpigz -p {threads} -c {input.ref} > {output.ref} "

#rule unzip_ref_gff_ncbi:
#    input:  ref = "{dir}/{species}/{ref}/annot/ncbi/{ref}.ncbi.gff.gz"
#    output: ref = "{dir}/{species}/{ref}/annot/ncbi/{ref}.ncbi.gff",
#    threads:  10
#    shell:
#        "unpigz -p {threads} -c {input.ref} > {output.ref} "


# rule copy_ref:
#   input:  IN_ANNOTATION_LIST
#   output: OUT_ANNOTATION_LIST
#   run:
#     for inp, out in zip(input, output):
#       shell(" ln {inp} {out} ")
