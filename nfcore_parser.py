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
    }
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
            if "enum" in details:
                param_entry["list"] = {item: item for item in details["enum"]}
                param_entry["type"] = "enum"
            if key == "input_output_options":
                output["gui_params"]["primary"][param] = param_entry
            else:
                output["gui_params"]["detailed"][param] = param_entry

# Save the output to the specified JSON file
with open(args.output_json, "w") as outfile:
    json.dump(output, outfile, indent=4)