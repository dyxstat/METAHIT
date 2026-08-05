import pandas as pd
import argparse

def merge_coverage(contig_info_path, coverage_path, output_path):
    coverage = pd.read_csv(coverage_path, sep="\t", header=0)
    if 'contigName' not in coverage.columns or 'totalAvgDepth' not in coverage.columns:
        raise ValueError("Coverage file must contain contigName and totalAvgDepth columns")

    coverage.rename(columns={'totalAvgDepth': 'coverage'}, inplace=True)
    coverage.rename(columns={'contigName': 'name'}, inplace=True)
    coverage = coverage[['name', 'coverage']]

    contig_info = pd.read_csv(contig_info_path)
    if 'name' not in contig_info.columns:
        raise ValueError("Contig info file must contain a name column")

    merged = contig_info.merge(coverage, on='name', how='left')
    merged.to_csv(output_path, index=False)
    print(f"[INFO] Merged contig info with coverage written to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--contig_info', required=True, help='Path to contig_info.csv')
    parser.add_argument('--coverage', required=True, help='Path to coverage.txt')
    parser.add_argument('--output', required=True, help='Path to output merged CSV')
    args = parser.parse_args()

    merge_coverage(args.contig_info, args.coverage, args.output)
