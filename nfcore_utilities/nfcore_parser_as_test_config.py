import json
import argparse
import os
import subprocess
import re

def parse_config_file(config_file):
    """
    Parse the Nextflow config file to extract parameters from the 'params' block.
    Returns a dictionary of parameter names and their default values.
    """
    params = {}
    inside_params_block = False

    try:
        with open(config_file, "r") as file:
            for line in file:
                # Detect the start of the 'params' block
                if line.strip().startswith("params {"):
                    inside_params_block = True
                    continue

                # Detect the end of the 'params' block
                if inside_params_block and line.strip() == "}":
                    inside_params_block = False
                    break

                # Extract parameters and their default values
                if inside_params_block:
                    match = re.match(r"(\w+)\s*=\s*(.+)", line.strip())
                    if match:
                        param_name = match.group(1)
                        param_value = match.group(2).split("//")[0].strip()  # Remove inline comments
                        params[param_name] = param_value
    except FileNotFoundError:
        print(f"Error: Config file '{config_file}' not found.")
    except Exception as e:
        print(f"Error while parsing config file: {e}")

    return params

   
# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert Nextflow schema JSON to a structured format.")
parser.add_argument("input_folder", help="Path to the input folder containing the required files.")
parser.add_argument("output_json", help="Path to the output JSON file.")
args = parser.parse_args()

# Extract workflow_type from the input folder name
workflow_type = os.path.basename(os.path.normpath(args.input_folder))


# Define the expected file paths
schema_file = os.path.join(args.input_folder, "nextflow_schema.json")
config_file = os.path.join(args.input_folder, "nextflow.config")
samplesheet_file = os.path.join(args.input_folder, "assets","samplesheet.csv")
readme_file = os.path.join(args.input_folder, "README.md")
usage_file = os.path.join(args.input_folder, "docs", "usage.md")

repo_path = args.input_folder  # Assuming the repository is in the input folder

# Define the output structure
output = {
    "task_name": "test",
    "entity_name": "test_nfcore",
    "globalResources": "/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/resources",
    "globalTmpdPath": "/tmp",
    "computing_type": "kubernetes",
    "Filesender_credential_file": "/mnt/data/ceitec_cfg2/710000-CEITEC/713000-cmm/713016-bioit/base/references_backup/resources_info/.secret/Filesender_credentials.json",
    "globalTaskPath": "/sequia",
    "S3_credential_file": "S3_credential_file",
    "aws_id": "aws_id",
    "aws_key": "aws_key"
}

# Check if the required files exist
# Parsing the main parameters:
if not os.path.isfile(schema_file):
    print(f"{schema_file} does not exist")
    print(f"Extract parameters from {config_file}")
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"{config_file} does not exist, cannot extract pipeline parameters")

    # Parse the config file
    config_params = parse_config_file(config_file)

    # Populate GUI parameters from config
    # output["gui_params"] = {"primary": {}, "detailed": {}}
    for param, value in config_params.items():
        if param == "input":
            output["input"] = "samplesheet.csv"
        elif param == "outdir":
            output["outdir"] = "results/" + workflow_type + "/"
        elif param == "genome":
            output["organism"] = "homo_sapiens"
            output["assembly"] = "GRCh38"
            output["release"] = "GRCh38_r111"
        else:
            output[param] = value.strip("'").strip('"')  # Remove quotes around strings

else:
    # Load the JSON schema
    print(f"Extract parameters from {schema_file}...")
    with open(schema_file, "r") as file:
        schema = json.load(file)

    # Populate GUI parameters from schema
    for key, value in (schema.get("$defs", {}) or schema.get("definitions", {})).items():
        if "properties" in value:
            for param, details in value["properties"].items():
                if param == "input":
                    output["input"] = "samplesheet.csv"
                elif param == "outdir":
                    output["outdir"] = "results/" + workflow_type + "/"
                elif param == "genome":
                    output["organism"] = "homo_sapiens"
                    output["assembly"] = "GRCh38"
                    output["release"] = "GRCh38_r111"
                else:
                    param_entry = {
                        "label": details.get("description", param),
                        "type": details.get("type", "string"),
                        "default": details.get("default", None),
                        "info": details.get("help_text", ""),
                    }
                    if param_entry["type"] == "string" and param_entry["default"] is None:
                        param_entry["default"] = ""
                    elif param_entry["type"] == "boolean" and param_entry["default"] is None:
                        param_entry["default"] = False

                    if "enum" in details:
                        param_entry["list"] = {item: item for item in details["enum"]}
                        param_entry["type"] = "enum"
                    if key == "input_output_options":
                        output[param] = param_entry["default"]
                    else:
                        output[param] = param_entry["default"]


def extract_csv_from_lines(lines):
    """
    Extract CSV headers from lines based on patterns.
    """
    for i, line in enumerate(lines):
        # Check for ```csv pattern
        if "```csv" in line:
            if i + 1 < len(lines):
                return lines[i + 1].strip().split(",")
        # Check for table header pattern (| delimiter)
        if line.strip().startswith("|") and line.strip().endswith("|"):
            # Split the line by '|' and remove empty entries
            return [col.strip() for col in line.strip().split("|") if col.strip()]
    return None

# Parsing the samplesheet structure:
csv_lines = []
found_pattern = False

# Check if samplesheet_file exists
if os.path.isfile(samplesheet_file):
    print(f"Parsing samplesheet structure from {samplesheet_file}...")
    with open(samplesheet_file, "r") as samplesheet:
        # Read the first line of the CSV file
        first_line = samplesheet.readline().strip()
        csv_lines = first_line.split(",")
        found_pattern = True
else:
    print(f"{samplesheet_file} does not exist")
    print(f"trying to parse samplesheet structure from {readme_file}...")
    if os.path.isfile(readme_file):
        # Search in README.md
        with open(readme_file, "r") as readme:
            lines = readme.readlines()
            csv_lines = extract_csv_from_lines(lines)
            if csv_lines:
                found_pattern = True
    else:
        print(f"{readme_file} does not exist")
        print(f"trying to parse samplesheet structure from {usage_file}...")

    # If not found in README.md, search in usage_file
    if not found_pattern and os.path.isfile(usage_file):
        with open(usage_file, "r") as usage:
            lines = usage.readlines()
            csv_lines = extract_csv_from_lines(lines)
            if csv_lines:
                found_pattern = True

# If still not found, print a message
if not found_pattern:
    print("No valid samplesheet structure found in samplesheet_file, README.md, or usage_file.")

# Process the extracted CSV lines
csv_json = {}
requested_params = []

if csv_lines:
    # Initialize the samples structure
    output["samples"] = {
        "111": {"sample_name": "first_sample"},
        "222": {"sample_name": "second_sample"},
        "333": {"sample_name": "third_sample"},
        "444": {"sample_name": "fourth_sample"},
    }

    for var in csv_lines:
        var = var.strip()
        if var in ["sample", "sample_id", "fastq_1", "filename_R1"]:
            continue
        elif var in ["paired", "fastq_2", "filename_R2"]:
            output["is_paired"] = True
        elif var == "strandedness":
            output["strandness"] = "reverse"
        else:
            # Add other elements with a default value of an empty string
            for sample_id in output["samples"]:
                output["samples"][sample_id][var] = ""

# Final output structure
final_output = {
    "task_name": output["task_name"],
    "entity_name": output["entity_name"],
}

# Add parsed parameters (excluding samples and global parameters)
for key, value in output.items():
    if key not in ["task_name", "entity_name", "samples", "globalResources", "globalTmpdPath", "computing_type",
                   "Filesender_credential_file", "globalTaskPath", "S3_credential_file", "aws_id", "aws_key"]:
        final_output[key] = value

# Add sample-related parameters
if "samples" in output:
    final_output["samples"] = output["samples"]

# Add global parameters
final_output.update({
    "globalResources": output["globalResources"],
    "globalTmpdPath": output["globalTmpdPath"],
    "computing_type": output["computing_type"],
    "Filesender_credential_file": output["Filesender_credential_file"],
    "globalTaskPath": output["globalTaskPath"],
    "S3_credential_file": output["S3_credential_file"],
    "aws_id": output["aws_id"],
    "aws_key": output["aws_key"],
})

# Save the final output to the specified JSON file
with open(args.output_json, "w") as outfile:
    json.dump(final_output, outfile, indent=4)