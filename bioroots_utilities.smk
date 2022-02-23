import os
import json
import pandas as pd


##### Reference processing #####
#
def global_ref_path():
  #"/mnt/references"
  return "/home/BioRoots/raw_fastq_qc/"


# setting reference
def load_ref(global_ref_path, config):
  if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(global_ref_path,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
  return config


# setting organism from reference
def load_organism(global_ref_path, config):
  f = open(os.path.join(global_ref_path,"reference_info","reference.json"),)
  reference_dict = json.load(f)
  f.close()
  config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
  return config


def reference_directory(global_ref_path, config):
  return os.path.join(global_ref_path,config["organism"],config["reference"])


##### Config processing #####
#
def load_sample(config):
  return pd.DataFrame.from_dict(config["samples"],orient="index")


def is_paired(config):
  if not config["is_paired"]:
    read_pair_tags = ["SE"]
    paired = "SE"
  else:
    read_pair_tags = ["R1", "R2"]
    paired = "PE"
  return [read_pair_tags, paired]
