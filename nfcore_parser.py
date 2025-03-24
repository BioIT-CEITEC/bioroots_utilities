import json

#filepath = "/Users/alejandromedaglia/Documents/nfcore_hackaton/rnaseq/nextflow_schema.json"
# Load the JSON schema
with open("nextflow_schema.json", "r") as file:
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

# Save the output to a new JSON file
with open("converted_schema.json", "w") as outfile:
    json.dump(output, outfile, indent=4)
