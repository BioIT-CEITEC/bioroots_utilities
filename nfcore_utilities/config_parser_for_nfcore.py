import os
import re
import json
import argparse

def find_file_with_structure(base_dir):
    # Define the target files to search for in priority order
    target_files = ["assets/samplesheet.csv", "README.md", "docs/usage.md"]
    found_files = []

    # Walk through the specified directory and subdirectories
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            relative_path = os.path.relpath(os.path.join(root, file), base_dir).replace("\\", "/")
            if relative_path in target_files:
                found_files.append(relative_path)  # Store relative paths

    # Return files in the order of priority
    return sorted(found_files, key=lambda x: target_files.index(x))

def extract_samplesheet_structure(file_path):
    structure = None
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()

            # Check if the file is a CSV file
            if file_path.endswith("samplesheet.csv"):
                structure = content.splitlines()[0]  # Assume the first line is the header

            # Check if the file is a README or usage file
            elif file_path.endswith(".md"):
                # Look for a line following ```csv
                match_csv_block = re.search(r"```csv\s*\n(.*?)(\n```|$)", content, re.DOTALL)
                if match_csv_block:
                    structure = match_csv_block.group(1).strip()
                else:
                    # Look for a line starting with |
                    match_table = re.search(r"^\|.*$", content, re.MULTILINE)
                    if match_table:
                        structure = match_table.group(0).strip()
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

    return structure

def parse_config_and_create_samplesheet(header, output_file):
    try:
        # Load the config.json file
        with open("config.json", "r", encoding="utf-8") as config_file:
            config_data = json.load(config_file)

        # Ensure the "samples" key exists in the config
        if "samples" not in config_data or not isinstance(config_data["samples"], dict):
            print("Error: 'samples' key not found or invalid in config.json.")
            return

        # Parse the header into columns
        columns = [col.strip() for col in header.split(",")]

        # Check if "sample" or "sample_id" is in the header
        sample_column = None
        if "sample" in columns:
            sample_column = "sample"
        elif "sample_id" in columns:
            sample_column = "sample_id"

        # Create the output samplesheet.csv file
        with open(output_file, "w", encoding="utf-8") as output_file:
            # Write the header
            output_file.write(header + "\n")

            # Write rows based on the config data
            for sample_id, sample_data in config_data["samples"].items():
                row = []
                for col in columns:
                    if col == sample_column:
                        # Use the sample_name value for "sample" or "sample_id" column
                        row.append(sample_data.get("sample_name", ""))
                    elif col in ["fastq_1", "filename_R1"]:
                        # Populate fastq_1 field
                        row.append(f"raw_fastq/{sample_data.get('sample_name', '')}_R1.fastq.gz")
                    elif col in ["fastq_2", "filename_R2"]:
                        # Populate fastq_2 field based on is_paired
                        if config_data.get("is_paired", False):
                            row.append(f"raw_fastq/{sample_data.get('sample_name', '')}_R2.fastq.gz")
                        else:
                            row.append("")
                    else:
                        # Use other column values or default to an empty string
                        row.append(str(sample_data.get(col, "")))
                output_file.write(",".join(row) + "\n")

        print(f"Output samplesheet created successfully at: {output_file.name}")

    except FileNotFoundError:
        print("Error: config.json file not found.")
    except json.JSONDecodeError:
        print("Error: Failed to parse config.json. Ensure it is valid JSON.")
    except Exception as e:
        print(f"Unexpected error: {e}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate a samplesheet.csv file from config.json.")
    parser.add_argument("folder", nargs="?", default="./", help="The folder to search for files (default: current directory).")
    parser.add_argument("output_file", nargs="?", default="samplesheet.csv", help="The output file name (default: samplesheet.csv).")
    args = parser.parse_args()

    base_dir = os.path.abspath(args.folder)
    output_file = args.output_file

    print(f"Searching for files in: {base_dir}")
    files = find_file_with_structure(base_dir)

    if not files:
        print("No target files found.")
        return

    # Prioritize extracting structure from the files in order
    for file_path in files:
        print(f"Found file: {file_path}")
        structure = extract_samplesheet_structure(file_path)
        if structure:
            print(f"Extracted structure from {file_path}:\n{structure}")
            parse_config_and_create_samplesheet(structure, output_file)
            return  # Stop after finding the first valid structure

    print("No valid structure found in any file.")

if __name__ == "__main__":
    main()