import json
import argparse
import os

# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert Nextflow schema JSON to a structured format.")
parser.add_argument("input_folder", help="Path to the input folder containing the required files.")
parser.add_argument("output_json", help="Path to the output JSON file.")
args = parser.parse_args()

# Define the expected file paths
schema_file = os.path.join(args.input_folder, "nextflow_schema.json")
readme_file = os.path.join(args.input_folder, "README.md")
usage_file = os.path.join(args.input_folder, "docs", "usage.md")

# Check if the required files exist
if not os.path.isfile(schema_file):
    raise FileNotFoundError(f"Required file not found: {schema_file}")
if not os.path.isfile(readme_file):
    raise FileNotFoundError(f"Required file not found: {readme_file}")
if not os.path.isfile(usage_file):
    raise FileNotFoundError(f"Required file not found: {usage_file}")

# Load the JSON schema
with open(schema_file, "r") as file:
    schema = json.load(file)

workflow_type = None
for key, value in schema.items():
    if key == "title" and isinstance(value, str) and value.startswith("nf-core/"):
        # Extract everything after "nf-core/" and before the next space
        workflow_type = value.split("nf-core/")[1].split()[0]
        break

# Define the output structure
output = {
    "workflow_description": {
        "name": workflow_type,
        "version": 1.0,
        "label": "workflow_type",
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

# Populate GUI parameters
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

# Open and parse the README file for ```csv patterns
csv_lines = []
found_pattern = False

with open(readme_file, "r") as readme:
    lines = readme.readlines()
    for i, line in enumerate(lines):
        if "```csv" in line:  # Check if the pattern is anywhere in the line
            found_pattern = True
            if i + 1 < len(lines):
                csv_lines.extend(lines[i + 1].strip().split(","))
            break

if not found_pattern:
    print("No ```csv pattern found in README.md. Searching in docs/usage.md...")
    with open(usage_file, "r") as usage:
        lines = usage.readlines()
        for i, line in enumerate(lines):
            if "```csv" in line:  # Check if the pattern is anywhere in the line
                found_pattern = True
                if i + 1 < len(lines):
                    csv_lines.extend(lines[i + 1].strip().split(","))
                break

if not found_pattern:
    print("No ```csv pattern found in either README.md or docs/usage.md.")

csv_json = {}
requested_params = []

for var in csv_lines:
    var = var.strip()
    if var in ["sample", "fastq_1", "fastq_2"]:
        continue
    if var == "paired":
        requested_params.append("is_paired")
    elif var == "strandedness":
        requested_params.append("strandness")
    else:
        csv_json[var] = {
            "label": var,
            "type": "string",
            "default": ""
        }

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