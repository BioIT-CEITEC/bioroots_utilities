#!/bin/bash

# Check if the correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_file> <output_folder>"
    exit 1
fi

# Assign arguments to variables
INPUT_FILE="$1"
OUTPUT_FOLDER="$2"

# Base URL for nf-core repositories
BASE_URL="https://github.com/nf-core"

# Check if the input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Create the output folder if it doesn't exist
if [[ ! -d "$OUTPUT_FOLDER" ]]; then
    mkdir -p "$OUTPUT_FOLDER"
fi

# Read the file line by line
while IFS= read -r repo_name || [[ -n "$repo_name" ]]; do
    # Skip empty lines
    if [[ -z "$repo_name" ]]; then
        continue
    fi

    # Construct the repository URL and directory path
    REPO_URL="$BASE_URL/$repo_name.git"
    REPO_DIR="$OUTPUT_FOLDER/$repo_name"

    # Check if the repository already exists
if [[ -d "$REPO_DIR" ]]; then
    echo "Repository '$repo_name' already exists. Checking for updates..."
    cd "$REPO_DIR"
    git fetch --all --tags

    # Determine the default branch (main or master)
    default_branch=$(git symbolic-ref refs/remotes/origin/HEAD | sed 's@^refs/remotes/origin/@@')

    latest_tag=$(git describe --tags `git rev-list --tags --max-count=1`)
    if [[ -n "$latest_tag" ]]; then
        echo "Checking out the latest release: $latest_tag"
        git checkout "$latest_tag"
    else
        echo "No tags found. Resetting to the latest commit on the default branch."
        git reset --hard "origin/$default_branch"
    fi
    cd ..
else
    # Clone the repository if it doesn't exist
    echo "Cloning repository: $repo_name"
    git clone "$REPO_URL" "$REPO_DIR"
    cd "$REPO_DIR"
    git fetch --all --tags

    # Determine the default branch (main or master)
    default_branch=$(git symbolic-ref refs/remotes/origin/HEAD | sed 's@^refs/remotes/origin/@@')

    latest_tag=$(git describe --tags `git rev-list --tags --max-count=1`)
    if [[ -n "$latest_tag" ]]; then
        echo "Checking out the latest release: $latest_tag"
        git checkout "$latest_tag"
    else
        echo "No tags found. Keeping the default branch."
        git checkout "$default_branch"
    fi
    cd ..
fi

    # Check if the operation was successful
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to process repository '$repo_name'."
    fi
done < "$INPUT_FILE"

echo "All repositories processed."