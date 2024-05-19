import os
import json
import pandas as pd
import boto3
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider



##### Check resource path
def check_resources():
    if "reference" in config:
        if "_r" in config["reference"]:
            globresource = "bioit"
            config["globalResources"] = config["globalResources"].replace("base/references_backup","resources")
        else:
            globresource = "bioda"
    else:
        if "references_backup" in config["globalResources"]:
            globresource = "bioda"
        else:
            globresource = "bioit"
    return globresource

##### Config processing #####
#
def load_sample():
    return pd.DataFrame.from_dict(config["samples"],orient="index")


def set_read_pair_tags():
  if not config["is_paired"]:
    read_pair_tags = [""]
  else:
    read_pair_tags = ["_R1", "_R2"]
  return read_pair_tags

def set_read_pair_qc_tags():
  if not config["is_paired"]:
    read_pair_qc_tags = ["SE"]
  else:
    read_pair_qc_tags = ["R1", "R2"]
  return read_pair_qc_tags

def set_paired_tags():
  if not config["is_paired"]:
    paired_tags = "SE"
  else:
    paired_tags = "PE"
  return paired_tags

def set_read_pair_dmtex_tags():
  if not config["is_paired"]:
    read_pair_dmtex_tags = ["_R1"]
  else:
    read_pair_dmtex_tags = ["_R1", "_R2"]
  return read_pair_dmtex_tags

##### kubernetes #####
##
#
if not "computing_type" in config:
  config["computing_type"] = "local"

if config["computing_type"] == "kubernetes":

    # with open(config["globalResources"] + "resources_info/.secret/S3_credentials.json") as f:
    #   S3_credentials = json.load(f)
    #S3 = S3RemoteProvider(host="https://storage-elixir1.cerit-sc.cz",access_key_id=S3_credentials["AWS_ID"],secret_access_key=S3_credentials["AWS_KEY"])
    S3 = S3RemoteProvider(host="https://storage-elixir1.cerit-sc.cz",access_key_id="acgt",secret_access_key="P84RsiL5TmHu0Ijd")
    client = boto3.client("s3",aws_access_key_id="acgt",aws_secret_access_key="P84RsiL5TmHu0Ijd",region_name="",endpoint_url="https://storage-elixir1.cerit-sc.cz")
    # S3_BUCKET = S3_credentials["S3_BUCKET"]
    S3_BUCKET = "acgt"
    task_directory = os.path.join(config["globalTaskPath"], config["task_name"]) + "/"

print(config["computing_type"])

##### Reference processing #####
##
#

####################
def load_dict(file_path):
    print(file_path)
    if config["computing_type"] == "kubernetes":
        if isinstance(file_path,list) and len(file_path) == 1:
            obj = client.get_object(Bucket=S3_BUCKET,Key=file_path[0])
            dictionary = json.loads(obj["Body"].read())
            return dictionary[0]
        else:
            if isinstance(file_path,str):
                obj = client.get_object(Bucket=S3_BUCKET,Key=file_path)
                dictionary = json.loads(obj["Body"].read())
                return dictionary
            else:
                obj = client.get_object(Bucket=S3_BUCKET,Key=file_path)
                dictionary = json.loads(obj["Body"].read())
                return (x for x in dictionary)
    else:
        if isinstance(file_path,list) and len(file_path) == 1:
            obj = open(file_path[0])
            dictionary = json.load(obj)
            obj.close()
            return dictionary[0]
        else:
            obj = open(file_path)
            dictionary = json.load(obj)
            obj.close()
            return dictionary


def load_ref():
    if config["lib_ROI"] != "wgs":
        # setting reference from lib_ROI
        lib_ROI_dict = load_dict(config["globalResources"] + "/reference_info/lib_ROI.json")
        config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
    return config


# def load_release():
#     release_cfg = load_dict(config["globalResources"] + "/reference_info/reference_release.json")
#     config["release"] = release_cfg[config["reference"]]["release"]
#     return config

def load_tooldir():
    globresource = check_resources()
    print(globresource)
    if globresource == "bioda":
        config["tooldir"] = os.path.join(config["globalResources"] ,"general")
    if globresource == "bioit":
        config["tooldir"] = os.path.join(config["globalResources"],"tools")


def load_ROI(globresource):
    k = open(os.path.join(config["globalResources"],"reference_info","lib_ROI.json"),)
    lib_ROI_dict = json.load(k)
    k.close()

    if config["lib_ROI"] == "rna":
        config["material"] = "RNA"
    else:
        if globresource == "bioda":
            config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if
                                   isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[
                                       ref_name].keys()][0]

        if globresource == "bioit":
            config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if
                                   isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[
                                       ref_name].keys()][0]
            config["lib_ROI"] = config["lib_ROI"].rsplit("_",1)[0]

    return config


def load_organism():
    globresource = check_resources()
    print(globresource)
    # setting organism from reference
    print(config["globalResources"])
    print(config["globalResources"] + "/reference_info/reference2.json")
    f = open(os.path.join(config["globalResources"],"reference_info","reference2.json"))
    reference_dict = json.load(f)
    f.close()
    k = open(os.path.join(config["globalResources"],"reference_info","kegg_reference.json"))
    kegg_dict = json.load(k)
    k.close()

    if "lib_ROI" in config and config["lib_ROI"] != "wgs":
        load_ROI(globresource)
 
    if globresource == "bioda":
        config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
        config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
        if len(config["species_name"].split(" (")) > 1:
            config["species"] = config["species_name"].split(" (")[1].replace(")","")
        config["reference_dir"] = os.path.join(config["globalResources"] , config["organism"] , config["reference"])
        config["organism_fasta"] = config["reference_dir"] + "/seq/" + config["reference"] + ".fa"
        config["organism_ucsc"] = config["reference_dir"] + "/seq/" + config["reference"] + ".fa.fai.ucsc"
        config["organism_gtf"] = config["reference_dir"] + "/annot/" + config["reference"] + ".gtf"
        config["organism_gtf_cellranger"] = config["reference_dir"] + "/annot/" + config["reference"] + "_cellranger.gtf"
        config["organism_cds_fasta"] = config["reference_dir"] + "/seq/" + config["reference"] + ".cds.fa"
        config["organism_cdna_fasta"] = config["reference_dir"] + "/seq/" + config["reference"] + ".cdna.fa"
        config["organism_star"] = config["reference_dir"] + "/index/STAR/SAindex"
        config["organism_starsolo"] = config["reference_dir"] + "/index/STAR/STAR_cellranger/SAindex"
        config["organism_rsem"] = config["reference_dir"] + "/index/RSEM/" + config["reference"] + ".idx.fa"
        config["organism_salmon"] = config["reference_dir"] + "/index/Salmon"
        config["organism_kallisto"] = config["reference_dir"] + "/index/Kallisto"
        config["organism_code"] = kegg_dict.get(config["species_name"])
        config["organism_picard_bed12"] = config["reference_dir"] + "/other/Picard_data/" + config["reference"] + ".bed12"
        config["organism_picard_refFlat"] = config["reference_dir"] + "/other/Picard_data/" + config["reference"] + ".refFlat"
        config["organism_ncbi_general"] = config["reference_dir"] + "/other/BOWTIE2/fastq_screen_RNA_indexes/" + config["reference"] + ".ncbi.fna",
        config["organism_ncbi_gff"] = config["reference_dir"] + "/other/BOWTIE2/fastq_screen_RNA_indexes/" + config["reference"] + ".ncbi.gff"
        config["organism_ncbi_rRNA"] = config["reference_dir"] + "/other/BOWTIE2/fastq_screen_RNA_indexes/" + config["reference"] + ".ncbi.rRNA.fasta"
        config["organism_ncbi_tRNA"] = config["reference_dir"] + "/other/BOWTIE2/fastq_screen_RNA_indexes/" + config["reference"] + ".ncbi.tRNA.fasta"
        config["organism_bwa"] = config["reference_dir"] + "/index/BWA/" + config["reference"] + ".bwt"
        config["organism_bowtie2"] = config["reference_dir"] + "/index/Bowtie2/" + config["reference"] + ".1.bt2"
        config["organism_vep_dir"] = config["reference_dir"] + "/annot/vep/"
        config["organism_chr_sizes"] = config["reference_dir"] + "/seq/" + config["reference"] + ".chrom.sizes"
        config["organism_dict"] = config["reference_dir"] + "/seq/" + config["reference"] + ".dict"
        config["organism_snp_bed"] = config["reference_dir"] + "/other/snp/" + config["reference"] + ".snp.bed"
        config["organism_custom_DB_folder"] = config["reference_dir"] + "/annot/custom_new2/"
        config["organism_cadd_db_snvs"] = config["reference_dir"] + "/annot/vep/CADD_scores_DB/whole_genome_SNVs.tsv.gz"
        config["organism_cadd_db_indels"] = config["reference_dir"] + "/annot/vep/CADD_scores_DB/gnomad.genomes.r3.0.indel.tsv.gz"
        config["organism_dbsnp"] = config["reference_dir"] + "/annot/dbSNP/common_all.vcf.gz"
        config["organism_cytoband"] = config["reference_dir"] + "/other/cytoband/" + config["reference"] + ".cytoband.tsv"        
        config["organism_svdb"] = config["reference_dir"] + "/other/svdb/gnomad_v2.1_sv.sites.vcf"
        config["organism_transcriptome"] = config["reference_dir"] + "/other/cellranger/refdata-gex-" + config["reference"]

    if globresource == "bioit":
        config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
        config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
        if len(config["species_name"].split(" (")) > 1:
            config["species"] = config["species_name"].split(" (")[1].replace(")","")
        config["assembly"] = config["reference"].rsplit("_",1)[0]
        config["release"] = config["reference"].rsplit("_",1)[1]
        config["reference_dir"] = os.path.join(config["globalResources"] , "references", config["organism"] , config["assembly"])
        config["organism_fasta"] = config["reference_dir"] + "/seq/" + config["assembly"] + ".fa"
        config["organism_ucsc"] = config["reference_dir"] + "/seq/" + config["assembly"] + ".fa.fai.ucsc"
        config["organism_gtf"] = config["reference_dir"] + "/annot/" + config["release"] + "/" + config["assembly"] + ".gtf"
        config["organism_gtf_cellranger"] = config["reference_dir"] + "/annot/" + config["release"] + "/" + config["assembly"] + "_cellranger.gtf"
        config["organism_cds_fasta"] = config["reference_dir"] + "/annot/" + config["release"] + "/" + config["assembly"] + ".cds.fa"
        config["organism_cdna_fasta"] = config["reference_dir"] + "/annot/" + config["release"] + "/" + config["assembly"] + ".cdna.fa"
        config["organism_star"] = config["reference_dir"] + "/tool_data/STAR/" + config["release"] + "/SAindex"
        config["organism_starsolo"] = config["reference_dir"] + "/tool_data/STAR/" + config["release"] + "/STAR_cellranger/SAindex"
        config["organism_rsem"] = config["reference_dir"] + "/tool_data/RSEM/" + config["release"] + "/" + config["assembly"] + ".idx.fa"
        config["organism_salmon"] = config["reference_dir"] + "/tool_data/Salmon/" + config["release"]
        config["organism_salmon_gentrome"] = config["reference_dir"] + "/tool_data/Salmon/" + config["release"] + "/Salmon_decoy/gentrome.fa"
        config["organism_kallisto"] = config["reference_dir"] + "/tool_data/Kallisto/" + config["release"] + "/Kallisto"
        config["organism_code"] = kegg_dict.get(config["species_name"])
        config["organism_picard_bed12"] = config["reference_dir"] + "/annot/" + config["release"] + "/Picard/" + config["assembly"] + ".bed12"
        config["organism_picard_refFlat"] = config["reference_dir"] + "/annot/" + config["release"] + "/Picard/" + config["assembly"] + ".refFlat"
        config["organism_ncbi_general"] = config["reference_dir"] + "/seq/BOWTIE2_fastq_screen/" + config["assembly"] + ".ncbi.fna",
        config["organism_ncbi_gff"] = config["reference_dir"] + "/seq/BOWTIE2_fastq_screen/" + config["assembly"] + ".ncbi.gff"
        config["organism_ncbi_rRNA"] = config["reference_dir"] + "/seq/BOWTIE2_fastq_screen/" + config["assembly"] + ".ncbi.rRNA.fasta"
        config["organism_ncbi_tRNA"] = config["reference_dir"] + "/seq/BOWTIE2_fastq_screen/" + config["assembly"] + ".ncbi.tRNA.fasta"
        config["organism_bwa"] = config["reference_dir"] + "/tool_data/BWA/" + config["assembly"] + ".bwt"
        config["organism_bowtie2"] = config["reference_dir"] + "/tool_data/Bowtie2/" + config["assembly"] + ".1.bt2"
        config["organism_vep_dir"] = config["reference_dir"] + "/annot/" + config["release"] + "/vep/"
        config["organism_chr_sizes"] = config["reference_dir"] + "/seq/" + config["assembly"] + ".chrom.sizes"
        config["organism_dict"] = config["reference_dir"] + "/seq/" + config["assembly"] + ".dict"
        config["organism_snp_bed"] = config["reference_dir"] + "/seq/" + config["assembly"] + ".snp.bed"
        config["organism_custom_DB_folder"] = config["reference_dir"] + "/annot/" + config["release"] + "/custom_new2/"
        config["organism_cadd_db_snvs"] = config["reference_dir"] + "/annot/" + config["release"] + "/vep/CADD_scores_DB/whole_genome_SNVs.tsv.gz"
        config["organism_cadd_db_indels"] = config["reference_dir"] + "/annot/" + config["release"] + "/vep/CADD_scores_DB/gnomad.genomes.r3.0.indel.tsv.gz"
        config["organism_dbsnp"] = config["reference_dir"] + "/annot/" + config["release"] + "/dbSNP/common_all.vcf.gz"
        config["organism_cytoband"] = config["reference_dir"] + "/others/cytoband/" + config["assembly"] + ".cytoband.tsv"
        config["organism_svdb"] = config["reference_dir"] + "/others/svdb/gnomad_v2.1_sv.sites.vcf"
        config["organism_transcriptome"] = config["reference_dir"] + "/tool_data/Cellranger/refdata-gex-" + config["assembly"]


    if "lib_ROI" in config:
        if globresource == "bioda":
            config["organism_dna_panel"] = config["reference_dir"] + "/intervals/" + config["lib_ROI"] + "/" + config[
                "lib_ROI"] + ".bed"
            config["organism_snps_panel"] = config["reference_dir"] + "/other/snp/" + config["lib_ROI"] + "/" + \
                                            config["lib_ROI"] + "_snps.bed"
            config["organism_interval_list"] = config["reference_dir"] + "/intervals/" + config["lib_ROI"] + "/" + \
                                               config["lib_ROI"] + ".interval_list"
        if globresource == "bioit":
            config["organism_dna_panel"] = config["reference_dir"] + "/others/DNA_ROI/" + config["lib_ROI"] + "/" + \
                                           config["lib_ROI"] + ".bed"
            config["organism_snps_panel"] = config["reference_dir"] + "/others/snp/" + config["lib_ROI"] + "/" + \
                                           config["lib_ROI"] + "_snps.bed"
            config["organism_interval_list"] = config["reference_dir"] + "/others/DNA_ROI/" + config["lib_ROI"] + "/" + \
                                               config["lib_ROI"] + ".interval_list"

    return config


def load_mirna():
    globresource = check_resources()
    print(globresource)
    # setting organism from reference
    print(config["globalResources"])
    print(config["globalResources"] + "/reference_info/reference2.json")
    f = open(os.path.join(config["globalResources"],"reference_info","reference2.json"),)
    reference_dict = json.load(f)
    f.close()
    k = open(os.path.join(config["globalResources"],"reference_info","kegg_reference.json"),)
    kegg_dict = json.load(k)
    k.close()

    if globresource == "bioda":
        config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
        config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
        if len(config["species_name"].split(" (")) > 1:
            config["species"] = config["species_name"].split(" (")[1].replace(")","")
        config["reference_dir"] = os.path.join(config["globalResources"] , config["organism"] , config["reference"])
        config["organism_rrna_star"] = config["reference_dir"] + "/index/STAR_rrna/SAindex"
        config["organism_mirbase"] = config["reference_dir"] + "/seq/hairpin.fa"
        config["organism_code"] = kegg_dict.get(config["species_name"])

    if globresource == "bioit":
        config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
        config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
        if len(config["species_name"].split(" (")) > 1:
            config["species"] = config["species_name"].split(" (")[1].replace(")","")
        config["assembly"] = config["reference"].rsplit("_",1)[1]
        config["reference_dir"] = os.path.join(config["globalResources"] , "references", config["organism"] , config["assembly"])
        config["organism_rrna_star"] = config["reference_dir"] + "/tool_data/STAR/SAindex"
        config["organism_mirbase"] = config["reference_dir"] + "/seq/hairpin.fa"
        config["organism_code"] = kegg_dict.get(config["species_name"])

    return config


def reference_directory():
    return os.path.join(config["globalResources"],config["organism"],config["reference"])


def remote_input_dir(dir_path: str):
    if isinstance(dir_path, list):
        directories = dir_path
    elif isinstance(dir_path, str):
        directories = [dir_path]
    contents = []
    if config["computing_type"] == "kubernetes":
        for path in directories:
            response = client.list_objects_v2(Bucket=S3_BUCKET, Prefix=path)
            contents += [S3.remote(os.path.join(S3_BUCKET, file_path["Key"])) for file_path in response["Contents"]]
    else:
        for path in directories:
            for root, dirs, files in os.walk(path,followlinks=True):
                for file in files:
                    contents.append(os.path.join(root,file))
    return contents


def get_path(filename):
    if len(filename) == 0:
        return filename
    if config["computing_type"] == "kubernetes":
        if os.path.isabs(filename[0]):
            if isinstance(filename,list) and len(filename) == 1:
                return S3_BUCKET + filename[0]
            else:
                if isinstance(filename,str):
                    return S3_BUCKET + filename
                return [S3_BUCKET + x for x in filename]

        else:
            if isinstance(filename,list) and len(filename) == 1:
                return S3_BUCKET + task_directory + filename[0]
            else:
                if isinstance(filename,str):
                    return S3_BUCKET + task_directory + filename
                return [S3_BUCKET + task_directory + x for x in filename]
    else:
        if isinstance(filename,list) and len(filename) == 1:
            return filename[0]
        return filename


def kubernetes_remote(remote_path):
     if isinstance(remote_path, list):
         return [S3.remote(file_path) for file_path in remote_path]
     return S3.remote(remote_path)


def remote(file_path):
    file_path = get_path(file_path)
    if config["computing_type"] == "kubernetes":
        return kubernetes_remote(file_path)
    return file_path


def get_bucket_name():
    if config["computing_type"] == "kubernetes":
        return S3_BUCKET
    return ""


##### Helper functions #####
##
# debugging function for listing all attributes and their classes in given snakemake object
def check_snakemake_object(snakemake, output_filename=None):
  common_atributes = ["append", "clear", "copy", "count", "extend", "get", "index", "insert", "items", "keys", "pop",
                      "remove", "reverse", "size", "size_mb", "sort"]

  if output_filename:
    original_stdout = sys.stdout  # Save a reference to the original standard output
    f = open(output_filename,"w")
    sys.stdout = f  # Change the standard output to the file we created.

  print("snakemake.inputs\n")
  for attr_name in [a for a in dir(snakemake.input) if not a.startswith("_") and not a in common_atributes]:
    print(attr_name)
    print(type(getattr(snakemake.input,attr_name)))
    print(getattr(snakemake.input,attr_name))
    print()
  print("-------------------------\n\nsnakemake.outputs\n")
  for attr_name in [a for a in dir(snakemake.output) if not a.startswith("_") and not a in common_atributes]:
    print(attr_name)
    print(type(getattr(snakemake.output,attr_name)))
    print(getattr(snakemake.output,attr_name))
    print()
  print("-------------------------\n\nsnakemake.params\n")
  for attr_name in [a for a in dir(snakemake.params) if not a.startswith("_") and not a in common_atributes]:
    print(attr_name)
    print(type(getattr(snakemake.params,attr_name)))
    print(getattr(snakemake.params,attr_name))
    print()

  if output_filename:
    sys.stdout = original_stdout  # Reset the standard output to its original value
    f.close()
