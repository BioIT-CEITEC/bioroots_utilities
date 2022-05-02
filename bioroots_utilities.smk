import os
import json
import pandas as pd
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

##### Reference processing #####
##
#


##### Config processing #####
#
def load_sample():
  return pd.DataFrame.from_dict(config["samples"],orient="index")


def set_read_pair_tags():
  if not config["is_paired"]:
    read_pair_tags = [""]
    paired = "SE"
  else:
    read_pair_tags = ["_R1", "_R2"]
    paired = "PE"
  return [read_pair_tags, paired]

##### kubernetes #####
##
#
if not "computing_type" in config:
  config["computing_type"] = "kubernetes"

if config["computing_type"] == "kubernetes":

  # with open(config["globalResources"] + "resources_info/.secret/S3_credentials.json") as f:
  #   S3_credentials = json.load(f)
  #S3 = S3RemoteProvider(host="https://storage-elixir1.cerit-sc.cz",access_key_id=S3_credentials["AWS_ID"],secret_access_key=S3_credentials["AWS_KEY"])
  S3 = S3RemoteProvider(host="https://storage-elixir1.cerit-sc.cz",access_key_id="acgt",secret_access_key="P84RsiL5TmHu0Ijd")
  # S3_BUCKET = S3_credentials["S3_BUCKET"]
  S3_BUCKET = "acgt"

####################
def load_ref():
  if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI

    f = open(os.path.join(config["biorootsResPath"],"resources_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
  return config


#setting organism from reference
def load_organism():
  f = open(os.path.join(config["biorootsResPath"],"resources_info","reference.json"))
  reference_dict = json.load(f)
  f.close()
  config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
  return config


def reference_directory():
  return os.path.join(config["globalResources"],"organisms",config["organism"],config["reference"])

####################


def remote(file_path):
  if config["computing_type"] == "kubernetes":
    #path = "/sequia/" + config["task_name"] + "/"

    if isinstance(file_path,list) and len(file_path) == 1:
      return S3.remote(S3_BUCKET + file_path[0])
    else:
      if isinstance(file_path,str):
        return S3.remote(S3_BUCKET + file_path)
      else:
        return S3.remote(S3_BUCKET + x for x in file_path)
  else:
    if isinstance(file_path,list) and len(file_path) == 1:
      return file_path[0]
    else:
      return file_path


##### Helper functions #####
##
# debugging function for listing all attributes and their classes in given snakemake object
def check_snakemake_object(snakemake, output_filename=None):
  common_atributes = ["append", "clear", "copy", "count", "extend", "get", "index", "insert", "items", "keys", "pop",
                      "remove", "reverse", "size", "size_mb", "sort"]

  if output_filename:
    original_stdout = sys.stdout  # Save a reference to the original standard output
    f = open(output_filename,'w')
    sys.stdout = f  # Change the standard output to the file we created.

  print("snakemake.inputs\n")
  for attr_name in [a for a in dir(snakemake.input) if not a.startswith('_') and not a in common_atributes]:
    print(attr_name)
    print(type(getattr(snakemake.input,attr_name)))
    print(getattr(snakemake.input,attr_name))
    print()
  print("-------------------------\n\nsnakemake.outputs\n")
  for attr_name in [a for a in dir(snakemake.output) if not a.startswith('_') and not a in common_atributes]:
    print(attr_name)
    print(type(getattr(snakemake.output,attr_name)))
    print(getattr(snakemake.output,attr_name))
    print()
  print("-------------------------\n\nsnakemake.params\n")
  for attr_name in [a for a in dir(snakemake.params) if not a.startswith('_') and not a in common_atributes]:
    print(attr_name)
    print(type(getattr(snakemake.params,attr_name)))
    print(getattr(snakemake.params,attr_name))
    print()

  if output_filename:
    sys.stdout = original_stdout  # Reset the standard output to its original value
    f.close()


rule all:
  input: load_ref()