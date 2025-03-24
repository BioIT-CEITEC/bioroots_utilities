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

# Define the output structure
output = {
    "workflow_description": {
        "name": "rnaseq_analysis",
        "version": 1.0,
        "label": "RNA-Seq Analysis",
        "type": "rnaseq_analysis",
        "inputs": "raw_fastq/{sample}*fastq.gz",
        "outputs": [
            "qc_reports/*",
            "logs/*",
            "sequences_summary/*",
            "aligned_reads/{sample}*aligned.bam"
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