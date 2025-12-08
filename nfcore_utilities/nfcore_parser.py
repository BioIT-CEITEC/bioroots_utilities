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

def get_git_version(repo_path):
    """
    Retrieve the version of the Git repository at the given path.
    Returns the latest tag or commit hash if no tags are available.
    """
    try:
        # Run 'git describe' to get the latest tag or commit hash
        result = subprocess.run(
            ["git", "-C", repo_path, "describe", "--tags", "--always"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error retrieving Git version: {e.stderr}")
        return None
    
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
git_version = get_git_version(repo_path)

# Define the output structure
output = {
    "workflow_description": {
        "name": workflow_type,
        "version": 1.0,
        "label": workflow_type,
        "type": "rnaseq_analysis",
        "inputs": "raw_fastq/{sample}*fastq.gz",
        "outputs": [
            "results/" + workflow_type + "/*"
        ],
        "report_index": "results/" + workflow_type + "/pipeline_info/pipeline_report.html",
        "reports": [
            "results/" + workflow_type + "/pipeline_info/pipeline_report.html"
        ]
    },
    "general_params": [
        "entity_name",
        "sample_name"
    ],
    "gui_params": {
        "primary": {},
        "detailed": {}
    },
}

if git_version:
    output["workflow_description"]["version"] = git_version
else:
    output["workflow_description"]["version"] = "1.0"

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
    output["gui_params"] = {"primary": {}, "detailed": {}}
    for param, value in config_params.items():
        if param == "input":
            output["gui_params"]["primary"]["input"] = {
                "type": "constant",
                "default": "samplesheet.csv"
            }
        elif param == "outdir":
            output["gui_params"]["primary"]["outdir"] = {
                "type": "constant",
                "default": "results/" + workflow_type + "/"
            }
        elif param == "genome":
            output["gui_params"]["primary"]["organism"] = {
                "label": "Organism",
                "type": "enum",
                "dynamicEnumName": "organism"
            }
            output["gui_params"]["primary"]["assembly"] = {
                "label": "Assembly",
                "type": "enum",
                "dynamicEnumName": "assembly",
                "filters": {
                    "group": {
                        "param": "organism",
                        "type": "value",
                        "showGroupLabel": False
                    }
                }
            }
            output["gui_params"]["primary"]["release"] = {
                "label": "Release",
                "type": "enum",
                "dynamicEnumName": "release",
                "filters": {
                    "group": {
                        "param": "assembly",
                        "type": "value",
                        "showGroupLabel": False
                    }
                }
            }
        else:
            output["gui_params"]["detailed"][param] = {
                "label": param,
                "type": "string" if value == "null" else "boolean" if value in ["true", "false"] else "string",
                "default": value.strip("'").strip('"')  # Remove quotes around strings
            }

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
                    output["gui_params"]["primary"]["input"] = {
                        "type": "constant",
                        "default": "samplesheet.csv"
                    }
                elif param == "outdir":
                    output["gui_params"]["primary"]["outdir"] = {
                        "type": "constant",
                        "default": "results/" + workflow_type + "/"
                    }
                elif param == "genome":
                    output["gui_params"]["primary"]["organism"] = {
                        "label": "Organism",
                        "type": "enum",
                        "dynamicEnumName": "organism"
                    }
                    output["gui_params"]["primary"]["assembly"] = {
                        "label": "Assembly",
                        "type": "enum",
                        "dynamicEnumName": "assembly",
                        "filters": {
                            "group": {
                                "param": "organism",
                                "type": "value",
                                "showGroupLabel": False
                            }
                        }
                    }
                    output["gui_params"]["primary"]["release"] = {
                        "label": "Release",
                        "type": "enum",
                        "dynamicEnumName": "release",
                        "filters": {
                            "group": {
                                "param": "assembly",
                                "type": "value",
                                "showGroupLabel": False
                            }
                        }
                    }
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
                        output["gui_params"]["primary"][param] = param_entry
                    else:
                        output["gui_params"]["detailed"][param] = param_entry


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
    for var in csv_lines:
        var = var.strip()
        if var in ["sample", "sample_id", "fastq_1", "filename_R1"]:
            continue
        if var in ["paired", "fastq_2", "filename_R2"]:
            requested_params.append("is_paired")
        elif var == "strandedness":
            requested_params.append("strandness")
        else:
            csv_json[var] = {
                "label": var,
                "type": "string",
                "default": ""
            }

# Final output structure
final_output = {
    "workflow_description": output["workflow_description"],
    "general_params": output["general_params"],
}

if requested_params:
    final_output["requested_params"] = requested_params

final_output["gui_params"] = output["gui_params"]

if csv_json:
    final_output["samples"] = csv_json

# Save the final output to the specified JSON file
with open(args.output_json, "w") as outfile:
    json.dump(final_output, outfile, indent=4)