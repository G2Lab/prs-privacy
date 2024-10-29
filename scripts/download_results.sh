#!/bin/bash

url="hhttps://github.com/G2Lab/prs-privacy-data/releases/download/results/results.zip"
destination_dir="results2"
archived_file="results.zip"

# Create destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Download the zip file
echo "Downloading $url..."
curl -L -o "$archived_file" "$url"  || { echo "Download failed"; exit 1; }

# Unzip the file to the destination directory
echo "Unpacking to $destination_dir..."
unzip -q -o "$archived_file" -d "$destination_dir" || { echo "Unzip failed"; exit 1; }
rm "$archived_file"

echo "Complete!"
