import pandas as pd


##### Config processing #####
def ref_path():
  return "/mnt/references"

def load_sample(config):
  sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
  return sample_tab

def is_paired(config):
  if not config["is_paired"]:
    read_pair_tags = ["SE"]
    paired = "SE"
  else:
    read_pair_tags = ["R1", "R2"]
    paired = "PE"
  return [read_pair_tags, paired]