import json
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Convert Nextflow schema JSON to a structured format.")
parser.add_argument("input_json", help="Path to the input JSON schema file.")
parser.add_argument("readme", help="Path to the README.md file (currently unused).")
parser.add_argument("output_json", help="Path to the output JSON file.")
args = parser.parse_args()

# Load the JSON schema
with open(args.input_json, "r") as file:
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
            "results/" + workflow_type+ "/*"
        ],
        "report_index": "qc_reports/multiqc/multiqc_report.html",
        "reports": [
            "qc_reports/multiqc/multiqc_report.html"
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
for key, value in schema.get("$defs", {}).items():
    if "properties" in value:
        for param, details in value["properties"].items():
            param_entry = {
                "label": details.get("description", param),
                "type": details.get("type", "string"),
                "default": details.get("default", None),
                "info": details.get("help_text", ""),
            }
            # Handle default values for specific types
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
with open(args.readme, "r") as readme_file:
    lines = readme_file.readlines()
    for i, line in enumerate(lines):
        if line.strip() == "```csv":
            # Capture the line following ```csv
            if i + 1 < len(lines):
                csv_lines.extend(lines[i + 1].strip().split(","))

# Convert csv_lines into a JSON structure with label, type, and default
csv_json = {
    var.strip(): {
        "label": var.strip(),
        "type": "string",
        "default": ""
    }
    for var in csv_lines
}

# Add the parsed CSV JSON structure to the output under "samples"
output["samples"] = csv_json

# Save the output to the specified JSON file
with open(args.output_json, "w") as outfile:
    json.dump(output, outfile, indent=4)