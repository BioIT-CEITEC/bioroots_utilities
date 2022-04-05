import os
import json
import pandas as pd
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

##### Reference processing #####
##
#

# setting reference
def load_ref():
  if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(config["globalResources"],"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
  return config


# setting organism from reference
def load_organism():
  f = open(os.path.join(config["globalResources"],"reference_info","reference.json"))
  reference_dict = json.load(f)
  f.close()
  config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
  return config


def reference_directory():
  return os.path.join(config["globalResources"],config["organism"],config["reference"])


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
if config["computing_type"] == "kubernetes":
  # f = open(config["globalResources"] + "/resources_info/S3_credentials.json")
  f = open(config["globalResources"] + "/reference_info/S3_credentials.json")
  S3_credentials = json.load(f)
  f.close()

  # S3 = S3RemoteProvider(host="https://storage-elixir1.cerit-sc.cz",access_key_id=S3_credentials["AWS_ID"],secret_access_key=S3_credentials["AWS_KEY"])
  # S3_BUCKET = S3_credentials["S3_BUCKET"]
  S3 = S3RemoteProvider(host="https://storage-elixir1.priv.cerit-sc.cz",access_key_id="acgt",secret_access_key="P84RsiL5TmHu0Ijd")
  S3_BUCKET = "acgt"

def remote(file_path):
  if config["computing_type"] == "kubernetes":
    path = "/sequia/" + config["task_name"] + "/"

    if isinstance(file_path,list) and len(file_path) == 1:
      return S3.remote(S3_BUCKET + path + file_path[0])
    else:
      if isinstance(file_path,str):
        return S3.remote(S3_BUCKET + path + file_path)
      else:
        return S3.remote(S3_BUCKET + path + x for x in file_path)
  else:
    if isinstance(file_path,list) and len(file_path) == 1:
      return file_path[0]
    else:
      return file_path

