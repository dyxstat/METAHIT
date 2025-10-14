import pandas as pd
import argparse

def merge_coverage(contig_info_path, coverage_path, output_path):
    # Read raw coverage file with all columns
    coverage = pd.read_csv(coverage_path, sep="\t", header=0)

    # Extract contig name from the composite field
    coverage['name'] = coverage['contigName'].str.extract(r'^(k\d+_\d+)')

    # Rename the coverage column
    coverage.rename(columns={'totalAvgDepth': 'coverage'}, inplace=True)

    # Keep only name and coverage
    coverage = coverage[['name', 'coverage']]

    # Read contig info
    contig_info = pd.read_csv(contig_info_path)

    # Merge by contig name
    merged = contig_info.merge(coverage, on='name', how='left')

    # Save
    merged.to_csv(output_path, index=False)
    print(f"[INFO] Merged contig info with coverage written to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--contig_info', required=True, help='Path to contig_info.csv')
    parser.add_argument('--coverage', required=True, help='Path to coverage.txt')
    parser.add_argument('--output', required=True, help='Path to output merged CSV')
    args = parser.parse_args()

    merge_coverage(args.contig_info, args.coverage, args.output)
