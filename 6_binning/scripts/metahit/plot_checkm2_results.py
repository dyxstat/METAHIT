#!/usr/bin/env python3
import sys
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

def main():
    # Check if correct number of arguments
    if len(sys.argv) != 3:
        print("Usage: python plot_checkm2_results.py <quality_report.tsv> <output_prefix>")
        sys.exit(1)
    
    quality_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    try:
        # Read CheckM2 quality report
        df = pd.read_csv(quality_file, sep='\t')
        
        # Sort by completeness (descending)
        df = df.sort_values('Completeness', ascending=False)
        
        # Extract bin names and metrics
        bins = df['Name'].tolist()
        completeness = df['Completeness'].tolist()
        contamination = df['Contamination'].tolist()
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, max(8, len(bins) * 0.3)))
        
        # Set width of bars
        width = 0.35
        
        # Set position of bars on x axis
        ind = np.arange(len(bins))
        
        # Create bars
        ax.barh(ind, completeness, width, color='#1f77b4', label='Completeness')
        ax.barh(ind + width, contamination, width, color='#ff7f0e', label='Contamination')
        
        # Add labels and title
        ax.set_xlabel('Percentage')
        ax.set_title('CheckM2 Quality Assessment')
        ax.set_yticks(ind + width/2)
        ax.set_yticklabels(bins)
        ax.legend()
        
        # Set limits
        ax.set_xlim(0, 105)
        
        # Add grid
        ax.grid(True, axis='x')
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f"{output_prefix}.png", dpi=300)
        
        # Create a second plot showing completeness vs contamination
        fig, ax = plt.subplots(figsize=(10, 8))
        scatter = ax.scatter(contamination, completeness, s=50, c=range(len(bins)), cmap='viridis')
        
        # Add labels for points
        for i, bin_name in enumerate(bins):
            ax.annotate(bin_name, (contamination[i], completeness[i]), 
                       xytext=(5, 0), textcoords='offset points', fontsize=8)
        
        # Add labels and title
        ax.set_xlabel('Contamination (%)')
        ax.set_ylabel('Completeness (%)')
        ax.set_title('CheckM2 Completeness vs Contamination')
        
        # Set limits
        ax.set_xlim(-1, max(contamination) + 5)
        ax.set_ylim(-1, 101)
        
        # Add grid
        ax.grid(True)
        
        # Save the scatter plot
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_scatter.png", dpi=300)
        
        print(f"Plots successfully generated: {output_prefix}.png and {output_prefix}_scatter.png")
        
    except Exception as e:
        print(f"Error generating plots: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()