#!/bin/bash
# Download the CSV file from a URL
wget "https://docs.google.com/spreadsheets/d/1uIritwhG4vXurJBDXyaERnZXObi75q5h3J5f309o-yQ/export?format=csv" -O "samplesheet.csv"
current_dir=$(pwd)

# Convert the CSV file to JSON using a Python script
python $current_dir/project/xsvato01/WES_CNV/CsvToJson.py $current_dir/samplesheet.csv $current_dir/samplesheet.json