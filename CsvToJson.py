import pandas as pd
import argparse
import json

# Create argument parser
parser = argparse.ArgumentParser(description='Convert CSV to JSON')
parser.add_argument('input_file', help='Input CSV filename')
parser.add_argument('output_file', help='Output JSON filename')
args = parser.parse_args()

# Read CSV and convert to JSON wrapped in a "runs" object
try:
    df = pd.read_csv(args.input_file)
    json_data = df.to_dict(orient='records')
    result = {"samples": json_data}
    with open(args.output_file, 'w') as json_file:
        json.dump(result, json_file, indent=4)
    print(f'Conversion successful. JSON file saved as {args.output_file}')
except Exception as e:
    print(f'Error: {e}')
