#!/usr/bin/env python3
"""
COMPREHENSIVE DENGUE SEVERITY IMMUNE REPERTOIRE ANALYSIS
===========================================================
This script performs comprehensive immune repertoire analysis for dengue severity
with integrated statistical visualizations. All results are saved in a RESULTS folder.

Features:
1. Diversity metrics calculation (richness, entropy, clonality)
2. Gini coefficient analysis for clonal expansion
3. Statistical comparisons between dengue severity groups
4. Publication-ready visualizations with statistical annotations
5. Multiple output formats (PDF, CSV, text reports)
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal, gaussian_kde
import itertools
import os
import re
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# ============================ CONFIGURATION ============================
# Create RESULTS directory
RESULTS_DIR = "RESULTS"
os.makedirs(RESULTS_DIR, exist_ok=True)

# Dengue severity metadata
DENGUE_METADATA = {
    'sample_id': [
        'gA1_S13', 'gA2_S14', 'gA3_S15', 'gA4_S16', 'gA5_S17', 'gA6_S18',
        'gB1_S7', 'gB2_S8', 'gB3_S9', 'gB4_S10', 'gB5_S11', 'gB6_S12',
        'gC1_S1', 'gC2_S2', 'gC3_S3', 'gC4_S4', 'gC5_S5', 'gC6_S6',
        'gD1_S19', 'gD2_S20', 'gD3_S21', 'gD4_S22', 'gD5_S23', 'gD6_S24'
    ],
    'group': [
        'Dengue_Shock_Syndrome', 'Dengue_Shock_Syndrome', 'Dengue_Shock_Syndrome',
        'Dengue_Shock_Syndrome', 'Dengue_Shock_Syndrome', 'Dengue_Shock_Syndrome',
        'Dengue_Hemorrhagic_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Hemorrhagic_Fever',
        'Dengue_Hemorrhagic_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Hemorrhagic_Fever',
        'Classical_Dengue_Fever', 'Classical_Dengue_Fever', 'Classical_Dengue_Fever',
        'Classical_Dengue_Fever', 'Classical_Dengue_Fever', 'Classical_Dengue_Fever',
        'Control', 'Control', 'Control', 'Control', 'Control', 'Control'
    ],
    'severity_level': [3]*6 + [2]*6 + [1]*6 + [0]*6
}

# Color palette for dengue severity
SEVERITY_PALETTE = {
    'Control': '#2E86AB',
    'Classical_Dengue_Fever': '#A23B72',
    'Dengue_Hemorrhagic_Fever': '#F18F01',
    'Dengue_Shock_Syndrome': '#C73E1D'
}

# ============================ PART 1: DIVERSITY METRICS ============================
print("="*70)
print("PART 1: DIVERSITY METRICS CALCULATION")
print("="*70)

def extract_total_reads_from_report(report_file):
    """Extract total sequencing reads from alignment report"""
    try:
        with open(report_file, 'r') as f:
            content = f.read()
            match = re.search(r'Total sequencing reads:\s*(\d+)', content)
            if match:
                return int(match.group(1))
    except Exception as e:
        print(f"Error reading {report_file}: {e}")
    return None

def calculate_entropy(read_counts):
    """Calculate Shannon entropy"""
    frequencies = read_counts / np.sum(read_counts)
    frequencies = frequencies[frequencies > 0]
    entropy = -np.sum(frequencies * np.log2(frequencies))
    return entropy

def calculate_richness(chain_reads, total_rnaseq_reads):
    """Calculate richness as defined in the Frontiers paper"""
    return chain_reads / (total_rnaseq_reads + chain_reads)

def process_sample(sample_id, results_dir):
    """Process a single sample and calculate all metrics"""
    
    report_file = os.path.join(results_dir, f"{sample_id}.align.report.txt")
    total_reads = extract_total_reads_from_report(report_file)
    
    if total_reads is None:
        print(f"Warning: Could not find total reads for {sample_id}")
        return None
    
    sample_results = []
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD']
    
    for chain in chains:
        chain_file = os.path.join(results_dir, f"{sample_id}.clones_{chain}.tsv")
        
        if not os.path.exists(chain_file):
            continue
            
        try:
            df = pd.read_csv(chain_file, sep='\t')
            
            if len(df) == 0:
                continue
                
            chain_total_reads = df['readCount'].sum()
            richness = calculate_richness(chain_total_reads, total_reads)
            entropy = calculate_entropy(df['readCount'].values)
            n_clonotypes = len(df)
            
            max_entropy = np.log2(n_clonotypes) if n_clonotypes > 0 else 0
            clonality = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
            
            sample_results.append({
                'sample_id': sample_id,
                'chain': chain,
                'richness': richness,
                'entropy': entropy,
                'n_clonotypes': n_clonotypes,
                'total_chain_reads': chain_total_reads,
                'clonality': clonality,
                'total_rnaseq_reads': total_reads
            })
            
        except Exception as e:
            print(f"Error processing {chain_file}: {e}")
            continue
    
    return sample_results

def calculate_all_diversity_metrics(results_dir="."):
    """Calculate diversity metrics for all samples"""
    
    all_files = os.listdir(results_dir)
    sample_ids = set()
    
    for file in all_files:
        if file.endswith('.align.report.txt'):
            sample_id = file.replace('.align.report.txt', '')
            sample_ids.add(sample_id)
    
    print(f"Found {len(sample_ids)} samples to process")
    
    all_results = []
    
    for sample_id in sorted(sample_ids):
        print(f"Processing {sample_id}...")
        sample_results = process_sample(sample_id, results_dir)
        
        if sample_results:
            all_results.extend(sample_results)
    
    if all_results:
        results_df = pd.DataFrame(all_results)
        results_df['log_richness'] = np.log10(results_df['richness'] + 1e-10)
        
        output_path = os.path.join(RESULTS_DIR, "all_samples_diversity_metrics.csv")
        results_df.to_csv(output_path, index=False)
        print(f"Diversity metrics saved to {output_path}")
        
        return results_df
    else:
        print("No diversity results to save")
        return None

# ============================ PART 2: INTEGRATE WITH METADATA ============================
print("\n" + "="*70)
print("PART 2: INTEGRATING WITH METADATA")
print("="*70)

def integrate_with_metadata(diversity_df):
    """Integrate diversity metrics with dengue metadata"""
    
    metadata_df = pd.DataFrame(DENGUE_METADATA)
    combined_df = pd.merge(diversity_df, metadata_df, on='sample_id', how='left')
    
    output_path = os.path.join(RESULTS_DIR, "dengue_repertoire_analysis_combined.csv")
    combined_df.to_csv(output_path, index=False)
    print(f"Combined data saved to {output_path}")
    
    print(f"\n=== MERGE VERIFICATION ===")
    print(f"Total samples in diversity data: {len(diversity_df['sample_id'].unique())}")
    print(f"Total samples after merge: {len(combined_df['sample_id'].unique())}")
    
    print(f"\n=== SAMPLE DISTRIBUTION BY GROUP ===")
    group_counts = combined_df[['sample_id', 'group']].drop_duplicates()['group'].value_counts()
    print(group_counts)
    
    return combined_df

# ============================ PART 3: STATISTICAL ANALYSIS ============================
print("\n" + "="*70)
print("PART 3: STATISTICAL ANALYSIS")
print("="*70)

class StatisticalResultsLogger:
    """Class to save all statistical results"""
    
    def __init__(self, filename="gini_statistical_analysis_report.txt"):
        self.filename = os.path.join(RESULTS_DIR, filename)
        self.results = []
        
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = f"""DENGUE SEVERITY IMMUNE REPERTOIRE ANALYSIS
===========================================================
Analysis Date: {timestamp}
Number of Samples: 24
Severity Groups: Control (n=6), Classical_Dengue_Fever (n=6), 
                 Dengue_Hemorrhagic_Fever (n=6), Dengue_Shock_Syndrome (n=6)
Analysis Type: Gini Coefficient and Diversity Metrics
===========================================================

"""
        self.results.append(header)
    
    def add_section(self, title):
        """Add a new section to results"""
        self.results.append(f"\n{'='*60}\n{title}\n{'='*60}\n")
    
    def add_result(self, text):
        """Add a single result line"""
        self.results.append(text + "\n")
    
    def add_table(self, title, df):
        """Add a dataframe as table"""
        self.add_section(title)
        self.results.append(df.to_string() + "\n")
    
    def save_to_file(self):
        """Save all results to file"""
        with open(self.filename, 'w') as f:
            f.writelines(self.results)
        print(f"Statistical report saved to {self.filename}")

def perform_group_comparisons(combined_df):
    """Perform statistical comparisons between dengue severity groups"""
    
    results = []
    chains = combined_df['chain'].unique()
    
    for chain in chains:
        chain_data = combined_df[combined_df['chain'] == chain]
        
        groups = ['Control', 'Classical_Dengue_Fever', 
                  'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
        group_data = [chain_data[chain_data['group'] == g]['richness'].values 
                     for g in groups if g in chain_data['group'].values]
        
        if len(group_data) > 1:
            h_stat, p_value = kruskal(*group_data)
            results.append({
                'chain': chain,
                'test': 'Kruskal-Wallis',
                'metric': 'richness',
                'statistic': h_stat,
                'p_value': p_value
            })
        
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                group1_data = chain_data[chain_data['group'] == groups[i]]['richness']
                group2_data = chain_data[chain_data['group'] == groups[j]]['richness']
                
                if len(group1_data) > 0 and len(group2_data) > 0:
                    u_stat, p_value = mannwhitneyu(group1_data, group2_data)
                    results.append({
                        'chain': chain,
                        'test': f'Mann-Whitney_{groups[i]}_vs_{groups[j]}',
                        'metric': 'richness',
                        'statistic': u_stat,
                        'p_value': p_value
                    })
    
    results_df = pd.DataFrame(results)
    
    # Save to file
    output_path = os.path.join(RESULTS_DIR, "statistical_analysis_results.csv")
    results_df.to_csv(output_path, index=False)
    print(f"Statistical analysis results saved to {output_path}")
    
    return results_df

# ============================ PART 4: GINI COEFFICIENT ANALYSIS ============================
print("\n" + "="*70)
print("PART 4: GINI COEFFICIENT ANALYSIS")
print("="*70)

def calculate_gini_coefficient(values):
    """Calculate Gini coefficient for a distribution"""
    if len(values) == 0:
        return np.nan
    
    values = np.array(values)
    values = values[values > 0]
    n = len(values)
    
    if n == 0:
        return np.nan
    
    sorted_values = np.sort(values)
    index = np.arange(1, n + 1)
    gini = (np.sum((2 * index - n - 1) * sorted_values)) / (n * np.sum(sorted_values))
    
    return abs(gini)

def calculate_clonal_expansion_gini(clone_data):
    """Calculate Gini coefficient for clonal expansion"""
    if len(clone_data) == 0:
        return np.nan
    
    vertex_sizes = clone_data['readCount'].values
    return calculate_gini_coefficient(vertex_sizes)

def process_sample_for_gini(sample_id, chain_file, chain_type):
    """Process a single sample for Gini coefficient"""
    try:
        df = pd.read_csv(chain_file, sep='\t')
        
        if len(df) == 0:
            return None
        
        expansion_gini = calculate_clonal_expansion_gini(df)
        
        return {
            'sample_id': sample_id,
            'chain': chain_type,
            'expansion_gini': expansion_gini,
            'n_clones': len(df),
            'total_reads': df['readCount'].sum()
        }
        
    except Exception as e:
        print(f"Error processing {chain_file}: {e}")
        return None

def calculate_all_gini_coefficients(results_dir="."):
    """Calculate Gini coefficients for all samples"""
    
    all_files = os.listdir(results_dir)
    sample_ids = set()
    
    for file in all_files:
        if file.endswith('.clones_IGH.tsv'):
            sample_id = file.replace('.clones_IGH.tsv', '')
            sample_ids.add(sample_id)
    
    print(f"Found {len(sample_ids)} samples to process")
    
    all_results = []
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    
    for sample_id in sorted(sample_ids):
        print(f"Processing {sample_id}...")
        
        for chain in chains:
            chain_file = os.path.join(results_dir, f"{sample_id}.clones_{chain}.tsv")
            
            if os.path.exists(chain_file):
                result = process_sample_for_gini(sample_id, chain_file, chain)
                if result:
                    all_results.append(result)
    
    if all_results:
        results_df = pd.DataFrame(all_results)
        
        # Add metadata
        metadata_df = pd.DataFrame(DENGUE_METADATA)
        results_df = pd.merge(results_df, metadata_df, on='sample_id', how='left')
        
        output_path = os.path.join(RESULTS_DIR, "dengue_gini_coefficient_analysis.csv")
        results_df.to_csv(output_path, index=False)
        print(f"Gini coefficient analysis saved to {output_path}")
        
        return results_df
    else:
        print("No Gini coefficient results")
        return None

# ============================ PART 5: STATISTICAL VISUALIZATION FUNCTIONS ============================
print("\n" + "="*70)
print("PART 5: STATISTICAL VISUALIZATION FUNCTIONS")
print("="*70)

def add_statistical_annotation(ax, x1, x2, y, h, color='black', linewidth=1.5):
    """Add statistical significance annotation to plots"""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=linewidth, c=color)
    ax.text((x1+x2)*0.5, y+h, "*", ha='center', va='bottom', color=color, 
            fontweight='bold', fontsize=12)

def perform_pairwise_tests(combined_df, metric, chain):
    """Perform pairwise Mann-Whitney tests between groups"""
    chain_data = combined_df[combined_df['chain'] == chain]
    groups = chain_data['group'].unique()
    
    results = []
    for group1, group2 in itertools.combinations(groups, 2):
        data1 = chain_data[chain_data['group'] == group1][metric].dropna()
        data2 = chain_data[chain_data['group'] == group2][metric].dropna()
        
        if len(data1) > 1 and len(data2) > 1:
            stat, p_value = mannwhitneyu(data1, data2)
            results.append({
                'group1': group1,
                'group2': group2,
                'p_value': p_value,
                'significant': p_value < 0.05
            })
    
    return results

def create_diversity_plots_with_stats_pdf(combined_df):
    """Create publication-ready diversity plots with statistical annotations in PDF format"""
    
    # Set style for publication
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (16, 12)
    plt.rcParams['pdf.fonttype'] = 42  # Ensure text is editable in PDF
    plt.rcParams['font.size'] = 12
    
    # Use the dengue severity palette from main analysis
    severity_palette = SEVERITY_PALETTE
    
    # Create the main figure
    fig = plt.figure(figsize=(20, 16))
    
    # 1. Richness by group and chain
    ax1 = plt.subplot(2, 2, 1)
    sns.boxplot(data=combined_df, x='chain', y='log_richness', hue='group', 
                palette=severity_palette, ax=ax1)
    plt.title('B-cell and T-cell Richness by Dengue Severity\n(Statistical Significance: *p < 0.05)', 
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45)
    plt.ylabel('Log Richness', fontsize=12)
    plt.xlabel('Chain Type', fontsize=12)
    plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add statistical annotations for richness
    chains = combined_df['chain'].unique()
    y_max = combined_df['log_richness'].max()
    for chain in chains:
        test_results = perform_pairwise_tests(combined_df, 'log_richness', chain)
        chain_idx = list(chains).index(chain)
        
        for i, result in enumerate(test_results):
            if result['significant']:
                group1_idx = list(combined_df['group'].unique()).index(result['group1'])
                group2_idx = list(combined_df['group'].unique()).index(result['group2'])
                
                # Adjust x positions for boxplot spacing
                x1 = chain_idx - 0.2 + group1_idx * 0.1
                x2 = chain_idx - 0.2 + group2_idx * 0.1
                y_pos = y_max + 0.1 + (i * 0.15)
                
                add_statistical_annotation(ax1, x1, x2, y_pos, 0.05)
    
    # 2. Entropy by group and chain
    ax2 = plt.subplot(2, 2, 2)
    sns.boxplot(data=combined_df, x='chain', y='entropy', hue='group', 
                palette=severity_palette, ax=ax2)
    plt.title('B-cell and T-cell Diversity by Dengue Severity\n(Statistical Significance: *p < 0.05)', 
              fontsize=14, fontweight='bold')
    plt.xticks(rotation=45)
    plt.ylabel('Shannon Entropy', fontsize=12)
    plt.xlabel('Chain Type', fontsize=12)
    
    # Add statistical annotations for entropy
    y_max_entropy = combined_df['entropy'].max()
    for chain in chains:
        test_results = perform_pairwise_tests(combined_df, 'entropy', chain)
        chain_idx = list(chains).index(chain)
        
        for i, result in enumerate(test_results):
            if result['significant']:
                group1_idx = list(combined_df['group'].unique()).index(result['group1'])
                group2_idx = list(combined_df['group'].unique()).index(result['group2'])
                
                x1 = chain_idx - 0.2 + group1_idx * 0.1
                x2 = chain_idx - 0.2 + group2_idx * 0.1
                y_pos = y_max_entropy + 0.5 + (i * 0.3)
                
                add_statistical_annotation(ax2, x1, x2, y_pos, 0.2)
    
    # 3. B-cell richness with detailed stats
    ax3 = plt.subplot(2, 2, 3)
    b_cells = combined_df[combined_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    sns.boxplot(data=b_cells, x='chain', y='log_richness', hue='group', 
                palette=severity_palette, ax=ax3)
    plt.title('B-cell Richness by Dengue Severity\nwith Statistical Annotations', 
              fontsize=14, fontweight='bold')
    plt.ylabel('Log Richness', fontsize=12)
    plt.xlabel('B-cell Chain Type', fontsize=12)
    
    # Add detailed annotations for B-cells
    b_chains = ['IGH', 'IGK', 'IGL']
    y_max_b = b_cells['log_richness'].max()
    
    for chain in b_chains:
        test_results = perform_pairwise_tests(b_cells, 'log_richness', chain)
        chain_idx = b_chains.index(chain)
        
        for i, result in enumerate(test_results):
            if result['significant']:
                group1_idx = list(b_cells['group'].unique()).index(result['group1'])
                group2_idx = list(b_cells['group'].unique()).index(result['group2'])
                
                x1 = chain_idx - 0.2 + group1_idx * 0.1
                x2 = chain_idx - 0.2 + group2_idx * 0.1
                y_pos = y_max_b + 0.1 + (i * 0.15)
                
                add_statistical_annotation(ax3, x1, x2, y_pos, 0.05)
                
                # Add p-value text for significant results
                p_text = f"p={result['p_value']:.3f}"
                ax3.text((x1+x2)*0.5, y_pos + 0.07, p_text, ha='center', va='bottom', 
                        fontsize=8, color='red', fontweight='bold')
    
    # 4. T-cell richness with detailed stats
    ax4 = plt.subplot(2, 2, 4)
    t_cells = combined_df[combined_df['chain'].isin(['TRA', 'TRB'])]
    if len(t_cells) > 0:
        sns.boxplot(data=t_cells, x='chain', y='log_richness', hue='group', 
                    palette=severity_palette, ax=ax4)
        plt.title('T-cell Richness by Dengue Severity\nwith Statistical Annotations', 
                  fontsize=14, fontweight='bold')
        plt.ylabel('Log Richness', fontsize=12)
        plt.xlabel('T-cell Chain Type', fontsize=12)
        
        # Add detailed annotations for T-cells
        t_chains = ['TRA', 'TRB']
        y_max_t = t_cells['log_richness'].max() if len(t_cells) > 0 else 0
        
        for chain in t_chains:
            chain_data = t_cells[t_cells['chain'] == chain]
            if len(chain_data) > 0:
                test_results = perform_pairwise_tests(t_cells, 'log_richness', chain)
                chain_idx = t_chains.index(chain)
                
                for i, result in enumerate(test_results):
                    if result['significant']:
                        group1_idx = list(t_cells['group'].unique()).index(result['group1'])
                        group2_idx = list(t_cells['group'].unique()).index(result['group2'])
                        
                        x1 = chain_idx - 0.2 + group1_idx * 0.1
                        x2 = chain_idx - 0.2 + group2_idx * 0.1
                        y_pos = y_max_t + 0.1 + (i * 0.15)
                        
                        add_statistical_annotation(ax4, x1, x2, y_pos, 0.05)
                        
                        # Add p-value text
                        p_text = f"p={result['p_value']:.3f}"
                        ax4.text((x1+x2)*0.5, y_pos + 0.07, p_text, ha='center', va='bottom', 
                                fontsize=8, color='red', fontweight='bold')
    
    plt.tight_layout()
    
    # Save to RESULTS directory
    output_path = os.path.join(RESULTS_DIR, "dengue_repertoire_diversity_with_stats.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Diversity plot with statistical annotations saved as PDF: {output_path}")

def create_enhanced_clonality_plot_pdf(combined_df):
    """Create enhanced clonality visualization with group information in PDF format"""
    
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (16, 8)
    plt.rcParams['pdf.fonttype'] = 42
    
    # Create a pivot table for clonality
    clonality_pivot = combined_df.pivot_table(
        index='sample_id', 
        columns='chain', 
        values='clonality'
    )
    
    # Add group information
    sample_groups = combined_df[['sample_id', 'group']].drop_duplicates().set_index('sample_id')
    clonality_pivot['group'] = sample_groups['group']
    
    # Define severity order
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    clonality_pivot['group'] = pd.Categorical(clonality_pivot['group'], categories=severity_order)
    clonality_pivot = clonality_pivot.sort_values('group')
    
    # Separate group column for coloring
    groups = clonality_pivot['group']
    clonality_data = clonality_pivot.drop('group', axis=1)
    
    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Heatmap
    im = ax1.imshow(clonality_data.T, aspect='auto', cmap='viridis', interpolation='nearest')
    ax1.set_yticks(range(len(clonality_data.columns)))
    ax1.set_yticklabels(clonality_data.columns)
    ax1.set_xlabel('Samples', fontweight='bold')
    ax1.set_title('Clonality Heatmap\n(Color by Chain Type)', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax1, label='Clonality')
    
    # Add group separators
    unique_groups = groups.unique()
    current_pos = 0
    for group in unique_groups:
        group_size = (groups == group).sum()
        ax1.axvline(x=current_pos + group_size - 0.5, color='red', linestyle='--', linewidth=2)
        ax1.text(current_pos + group_size/2, -1, group, 
                ha='center', va='top', fontweight='bold', rotation=45)
        current_pos += group_size
    
    # Group-wise clonality summary
    group_clonality = combined_df.groupby(['group', 'chain'])['clonality'].mean().unstack()
    group_clonality = group_clonality.reindex(severity_order)
    sns.heatmap(group_clonality, annot=True, fmt='.3f', cmap='viridis', ax=ax2, 
                cbar_kws={'label': 'Mean Clonality'})
    ax2.set_title('Mean Clonality by Group and Chain', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Dengue Severity Group', fontweight='bold')
    
    plt.tight_layout()
    
    # Save to RESULTS directory
    output_path = os.path.join(RESULTS_DIR, "enhanced_clonality_analysis.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Enhanced clonality plot saved as PDF: {output_path}")

def create_individual_chain_plots_pdf_comprehensive(combined_df):
    """Create individual plots for each chain with detailed statistics in PDF format"""
    
    sns.set_style("whitegrid")
    plt.rcParams['pdf.fonttype'] = 42
    
    chains = combined_df['chain'].unique()
    n_chains = len(chains)
    
    fig, axes = plt.subplots(2, (n_chains + 1) // 2, figsize=(20, 12))
    axes = axes.flatten()
    
    severity_palette = SEVERITY_PALETTE
    
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    
    for idx, chain in enumerate(chains):
        if idx < len(axes):
            chain_data = combined_df[combined_df['chain'] == chain]
            
            # Richness plot
            sns.boxplot(data=chain_data, x='group', y='log_richness', 
                       palette=severity_palette, order=severity_order, ax=axes[idx])
            sns.stripplot(data=chain_data, x='group', y='log_richness', 
                         color='black', alpha=0.6, order=severity_order, ax=axes[idx])
            
            # Add statistical annotations
            test_results = perform_pairwise_tests(combined_df, 'log_richness', chain)
            y_max = chain_data['log_richness'].max()
            
            significant_pairs = []
            for i, result in enumerate(test_results):
                if result['significant']:
                    group1_idx = severity_order.index(result['group1'])
                    group2_idx = severity_order.index(result['group2'])
                    
                    x1, x2 = group1_idx, group2_idx
                    y_pos = y_max + 0.1 + (len(significant_pairs) * 0.15)
                    
                    add_statistical_annotation(axes[idx], x1, x2, y_pos, 0.05)
                    significant_pairs.append(f"{result['group1']} vs {result['group2']}: p={result['p_value']:.3f}")
            
            # Count number of significant pairs
            n_sig = len(significant_pairs)
            
            axes[idx].set_title(f'{chain} Chain - Richness\n{n_sig} significant comparison(s)', 
                               fontsize=12, fontweight='bold')
            axes[idx].set_ylabel('Log Richness')
            axes[idx].set_xlabel('')
            axes[idx].tick_params(axis='x', rotation=45)
            
            # Add p-values as text
            if significant_pairs:
                axes[idx].text(0.02, 0.98, '\n'.join(significant_pairs), 
                              transform=axes[idx].transAxes, va='top', 
                              fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Remove empty subplots
    for idx in range(len(chains), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.tight_layout()
    
    # Save to RESULTS directory
    output_path = os.path.join(RESULTS_DIR, "individual_chain_analysis.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Individual chain plots saved as PDF: {output_path}")

def create_statistical_summary_csv(combined_df):
    """Create a summary table of statistical results"""
    
    summary_data = []
    
    for chain in combined_df['chain'].unique():
        for metric in ['log_richness', 'entropy']:
            test_results = perform_pairwise_tests(combined_df, metric, chain)
            
            for result in test_results:
                if result['significant']:
                    summary_data.append({
                        'Chain': chain,
                        'Metric': metric,
                        'Comparison': f"{result['group1']} vs {result['group2']}",
                        'P-value': result['p_value'],
                        'Significance': 'Yes' if result['significant'] else 'No'
                    })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values(['Chain', 'Metric', 'P-value'])
        
        output_path = os.path.join(RESULTS_DIR, "statistical_significance_summary.csv")
        summary_df.to_csv(output_path, index=False)
        
        print("\n=== STATISTICAL SIGNIFICANCE SUMMARY ===")
        print(summary_df.to_string(index=False))
        print(f"\nSummary saved to: {output_path}")
    else:
        print("No statistically significant differences found.")

# ============================ PART 6: ADDITIONAL VISUALIZATIONS ============================
def create_diversity_plots_pdf(combined_df):
    """Create diversity plots in PDF format"""
    
    sns.set_style("whitegrid")
    plt.rcParams['pdf.fonttype'] = 42  # Ensure text is editable in PDF
    
    # 1. Main diversity plot
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Richness by group and chain
    ax1 = axes[0, 0]
    sns.boxplot(data=combined_df, x='chain', y='log_richness', hue='group', 
                palette=SEVERITY_PALETTE, ax=ax1)
    ax1.set_title('B-cell and T-cell Richness by Dengue Severity', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Log Richness')
    ax1.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Entropy by group and chain
    ax2 = axes[0, 1]
    sns.boxplot(data=combined_df, x='chain', y='entropy', hue='group', 
                palette=SEVERITY_PALETTE, ax=ax2)
    ax2.set_title('B-cell and T-cell Diversity by Dengue Severity', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Shannon Entropy')
    
    # B-cell richness
    ax3 = axes[1, 0]
    b_cells = combined_df[combined_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    sns.boxplot(data=b_cells, x='chain', y='log_richness', hue='group', 
                palette=SEVERITY_PALETTE, ax=ax3)
    ax3.set_title('B-cell Richness by Dengue Severity', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Log Richness')
    
    # T-cell richness
    ax4 = axes[1, 1]
    t_cells = combined_df[combined_df['chain'].isin(['TRA', 'TRB'])]
    if len(t_cells) > 0:
        sns.boxplot(data=t_cells, x='chain', y='log_richness', hue='group', 
                    palette=SEVERITY_PALETTE, ax=ax4)
        ax4.set_title('T-cell Richness by Dengue Severity', fontsize=14, fontweight='bold')
        ax4.set_ylabel('Log Richness')
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "dengue_repertoire_diversity.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Diversity plot saved as PDF: {output_path}")

def create_gini_plots_pdf(gini_df):
    """Create Gini coefficient plots in PDF format"""
    
    sns.set_style("whitegrid")
    
    # Main Gini analysis figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    
    # 1. Expansion Gini by chain and severity
    ax1 = axes[0, 0]
    sns.boxplot(data=gini_df, x='chain', y='expansion_gini', hue='group', 
                ax=ax1, palette=SEVERITY_PALETTE, order=['IGH', 'IGK', 'IGL', 'TRA', 'TRB'])
    ax1.set_title('Clonal Expansion in Dengue Severity', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Gini Coefficient')
    ax1.legend(title='Dengue Severity', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 2. B-cell chains only
    ax2 = axes[0, 1]
    b_cells = gini_df[gini_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    sns.boxplot(data=b_cells, x='chain', y='expansion_gini', hue='group', 
                ax=ax2, palette=SEVERITY_PALETTE, order=['IGH', 'IGK', 'IGL'])
    ax2.set_title('B-cell Clonal Expansion', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Gini Coefficient')
    
    # 3. T-cell chains only
    ax3 = axes[0, 2]
    t_cells = gini_df[gini_df['chain'].isin(['TRA', 'TRB'])]
    if len(t_cells) > 0:
        sns.boxplot(data=t_cells, x='chain', y='expansion_gini', hue='group', 
                    ax=ax3, palette=SEVERITY_PALETTE, order=['TRA', 'TRB'])
        ax3.set_title('T-cell Clonal Expansion', fontsize=14, fontweight='bold')
        ax3.set_ylabel('Gini Coefficient')
    
    # 4. Heatmap
    ax4 = axes[1, 0]
    gini_pivot = gini_df.pivot_table(index='sample_id', columns='chain', 
                                     values='expansion_gini')
    gini_pivot['group'] = gini_df[['sample_id', 'group']].drop_duplicates().set_index('sample_id')['group']
    gini_pivot = gini_pivot.sort_values('group')
    groups = gini_pivot['group']
    gini_data = gini_pivot.drop('group', axis=1)
    
    im = ax4.imshow(gini_data.T, aspect='auto', cmap='RdYlBu_r', 
                   vmin=0, vmax=1, interpolation='nearest')
    ax4.set_title('Clonal Expansion Gini Heatmap', fontsize=14, fontweight='bold')
    ax4.set_xlabel('Samples (Grouped by Severity)')
    ax4.set_ylabel('Immune Receptor Chains')
    plt.colorbar(im, ax=ax4, label='Gini Coefficient')
    
    # 5. Gini vs Severity Level
    ax5 = axes[1, 1]
    for chain in gini_df['chain'].unique():
        chain_data = gini_df[gini_df['chain'] == chain]
        ax5.scatter(chain_data['severity_level'], chain_data['expansion_gini'], 
                   alpha=0.6, label=chain, s=50, edgecolors='black')
    ax5.set_title('Gini Coefficient vs Dengue Severity', fontsize=14, fontweight='bold')
    ax5.set_xlabel('Dengue Severity Level')
    ax5.set_ylabel('Gini Coefficient')
    ax5.set_xticks([0, 1, 2, 3])
    ax5.set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'])
    ax5.legend(title='Chain Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 6. Distribution by severity group
    ax6 = axes[1, 2]
    sns.violinplot(data=gini_df, x='group', y='expansion_gini', 
                   ax=ax6, palette=SEVERITY_PALETTE, order=severity_order)
    ax6.set_title('Distribution of Gini Coefficients', fontsize=14, fontweight='bold')
    ax6.set_ylabel('Gini Coefficient')
    ax6.set_xlabel('Dengue Severity Group')
    ax6.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "dengue_gini_coefficient_analysis.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Gini analysis plot saved as PDF: {output_path}")

def create_bcell_tcell_comparison_pdf(gini_df):
    """Create B-cell vs T-cell comparison plot in PDF"""
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Separate B and T cells
    b_cells = gini_df[gini_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    t_cells = gini_df[gini_df['chain'].isin(['TRA', 'TRB'])]
    
    # Calculate mean Gini per group
    b_mean = b_cells.groupby('group')['expansion_gini'].mean()
    t_mean = t_cells.groupby('group')['expansion_gini'].mean()
    
    # Bar plot comparison
    x = np.arange(len(b_mean))
    width = 0.35
    
    axes[0].bar(x - width/2, b_mean.values, width, label='B-cells', 
               color='#3498db', alpha=0.8)
    axes[0].bar(x + width/2, t_mean.values, width, label='T-cells', 
               color='#e74c3c', alpha=0.8)
    
    axes[0].set_title('B-cell vs T-cell Clonal Expansion', fontsize=14, fontweight='bold')
    axes[0].set_ylabel('Mean Gini Coefficient')
    axes[0].set_xlabel('Dengue Severity Group')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(b_mean.index, rotation=45)
    axes[0].legend()
    
    # Percentage change from control
    axes[1].bar(x, ((b_mean - b_mean.loc['Control']) / b_mean.loc['Control'] * 100).values,
               color=[SEVERITY_PALETTE[g] for g in b_mean.index], alpha=0.8)
    
    axes[1].set_title('Clonal Expansion Change Relative to Control', fontsize=14, fontweight='bold')
    axes[1].set_ylabel('% Change from Control')
    axes[1].set_xlabel('Dengue Severity Group')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(b_mean.index, rotation=45)
    axes[1].axhline(y=0, color='black', linestyle='-', linewidth=1)
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "dengue_bcell_tcell_comparison.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"B-cell/T-cell comparison saved as PDF: {output_path}")

def create_enhanced_clonality_pdf(combined_df):
    """Create enhanced clonality plot in PDF"""
    
    clonality_pivot = combined_df.pivot_table(
        index='sample_id', 
        columns='chain', 
        values='clonality'
    )
    
    sample_groups = combined_df[['sample_id', 'group']].drop_duplicates().set_index('sample_id')
    clonality_pivot['group'] = sample_groups['group']
    clonality_pivot = clonality_pivot.sort_values('group')
    
    groups = clonality_pivot['group']
    clonality_data = clonality_pivot.drop('group', axis=1)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Heatmap
    im = ax1.imshow(clonality_data.T, aspect='auto', cmap='viridis', interpolation='nearest')
    ax1.set_title('Clonality Heatmap', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Samples')
    ax1.set_ylabel('Chain Type')
    plt.colorbar(im, ax=ax1, label='Clonality')
    
    # Group-wise clonality summary
    group_clonality = combined_df.groupby(['group', 'chain'])['clonality'].mean().unstack()
    sns.heatmap(group_clonality, annot=True, cmap='viridis', ax=ax2, 
                cbar_kws={'label': 'Mean Clonality'})
    ax2.set_title('Mean Clonality by Group and Chain', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Dengue Severity Group')
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "enhanced_clonality_analysis_comprehensive.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Enhanced clonality plot saved as PDF: {output_path}")

def create_individual_chain_plots_pdf(gini_df):
    """Create individual chain plots in PDF"""
    
    chains = gini_df['chain'].unique()
    n_chains = len(chains)
    
    fig, axes = plt.subplots(2, (n_chains + 1) // 2, figsize=(20, 12))
    axes = axes.flatten()
    
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    
    for idx, chain in enumerate(chains):
        if idx < len(axes):
            chain_data = gini_df[gini_df['chain'] == chain]
            
            # Boxplot with individual points
            sns.boxplot(data=chain_data, x='group', y='expansion_gini', 
                       palette=SEVERITY_PALETTE, order=severity_order, ax=axes[idx])
            sns.stripplot(data=chain_data, x='group', y='expansion_gini', 
                         color='black', alpha=0.6, order=severity_order, ax=axes[idx])
            
            axes[idx].set_title(f'{chain} Chain - Clonal Expansion', fontsize=12, fontweight='bold')
            axes[idx].set_ylabel('Gini Coefficient')
            axes[idx].set_xlabel('')
            axes[idx].tick_params(axis='x', rotation=45)
    
    # Remove empty subplots
    for idx in range(len(chains), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "individual_chain_analysis_gini.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Individual chain plots (Gini) saved as PDF: {output_path}")

def create_statistical_summary_pdf(gini_df):
    """Create statistical summary in PDF"""
    
    fig = plt.figure(figsize=(20, 15))
    
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    chain_order = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    
    # ========== SUBPLOT 1: All chains scatter ==========
    ax1 = plt.subplot(3, 3, 1)
    for chain in chain_order:
        chain_data = gini_df[gini_df['chain'] == chain]
        if len(chain_data) > 0:
            for i, group in enumerate(severity_order):
                group_data = chain_data[chain_data['group'] == group]
                if len(group_data) > 0:
                    jitter = np.random.normal(i, 0.1, len(group_data))
                    ax1.scatter(jitter, group_data['expansion_gini'], 
                              color=SEVERITY_PALETTE[group], alpha=0.6, s=50, 
                              label=group if chain == 'IGH' else "")
    
    ax1.set_title('All Chains: Clonal Expansion by Severity', fontweight='bold')
    ax1.set_xlabel('Dengue Severity Groups')
    ax1.set_ylabel('Gini Coefficient')
    ax1.set_xticks(range(len(severity_order)))
    ax1.set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'], rotation=45)
    
    # ========== SUBPLOT 2: B-cell statistical boxplot ==========
    ax2 = plt.subplot(3, 3, 2)
    b_cells = gini_df[gini_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    if len(b_cells) > 0:
        sns.boxplot(data=b_cells, x='group', y='expansion_gini', hue='group',
                   ax=ax2, palette=SEVERITY_PALETTE, order=severity_order, legend=False)
        sns.stripplot(data=b_cells, x='group', y='expansion_gini', 
                     ax=ax2, color='black', alpha=0.5, size=4, order=severity_order)
        
        ax2.set_title('B-cell Chains: Statistical Comparison', fontweight='bold')
        ax2.set_xlabel('Dengue Severity Groups')
        ax2.set_ylabel('Gini Coefficient')
        ax2.set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'], rotation=45)
    
    # ========== SUBPLOT 3: T-cell statistical boxplot ==========
    ax3 = plt.subplot(3, 3, 3)
    t_cells = gini_df[gini_df['chain'].isin(['TRA', 'TRB'])]
    if len(t_cells) > 0:
        sns.boxplot(data=t_cells, x='group', y='expansion_gini', hue='group',
                   ax=ax3, palette=SEVERITY_PALETTE, order=severity_order, legend=False)
        sns.stripplot(data=t_cells, x='group', y='expansion_gini', 
                     ax=ax3, color='black', alpha=0.5, size=4, order=severity_order)
        
        ax3.set_title('T-cell Chains: Statistical Comparison', fontweight='bold')
        ax3.set_xlabel('Dengue Severity Groups')
        ax3.set_ylabel('Gini Coefficient')
        ax3.set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'], rotation=45)
    
    # ========== SUBPLOT 4: Statistical significance heatmap ==========
    ax4 = plt.subplot(3, 3, 4)
    
    p_value_matrix = pd.DataFrame(index=chain_order, 
                                 columns=['CvsD', 'CvsH', 'CvsS', 'DvsH', 'DvsS', 'HvsS'])
    
    for chain in chain_order:
        chain_data = gini_df[gini_df['chain'] == chain]
        
        # Control vs Classical
        c_data = chain_data[chain_data['group'] == 'Control']['expansion_gini']
        d_data = chain_data[chain_data['group'] == 'Classical_Dengue_Fever']['expansion_gini']
        if len(c_data) > 1 and len(d_data) > 1:
            _, p = mannwhitneyu(c_data, d_data)
            p_value_matrix.loc[chain, 'CvsD'] = p
        
        # Control vs Hemorrhagic
        h_data = chain_data[chain_data['group'] == 'Dengue_Hemorrhagic_Fever']['expansion_gini']
        if len(c_data) > 1 and len(h_data) > 1:
            _, p = mannwhitneyu(c_data, h_data)
            p_value_matrix.loc[chain, 'CvsH'] = p
        
        # Control vs Shock
        s_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['expansion_gini']
        if len(c_data) > 1 and len(s_data) > 1:
            _, p = mannwhitneyu(c_data, s_data)
            p_value_matrix.loc[chain, 'CvsS'] = p
    
    p_values = p_value_matrix.astype(float).fillna(1)
    im = ax4.imshow(-np.log10(p_values), cmap='YlOrRd', aspect='auto', vmin=0, vmax=3)
    
    ax4.set_title('Statistical Significance Heatmap', fontweight='bold')
    ax4.set_xlabel('Group Comparisons')
    ax4.set_ylabel('Immune Receptor Chain')
    ax4.set_xticks(range(len(p_value_matrix.columns)))
    ax4.set_xticklabels(['CvsD', 'CvsH', 'CvsS', 'DvsH', 'DvsS', 'HvsS'], rotation=45)
    ax4.set_yticks(range(len(chain_order)))
    ax4.set_yticklabels(chain_order)
    plt.colorbar(im, ax=ax4, label='-log10(p-value)')
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "dengue_gini_statistical_summary.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Statistical summary saved as PDF: {output_path}")

# ============================ PART 7: MAIN EXECUTION ============================
print("\n" + "="*70)
print("MAIN ANALYSIS PIPELINE")
print("="*70)

def main():
    """Main analysis pipeline"""
    
    print("Starting comprehensive immune repertoire analysis...")
    
    # Part 1: Calculate diversity metrics
    print("\n1. Calculating diversity metrics...")
    diversity_df = calculate_all_diversity_metrics()
    
    if diversity_df is None:
        print("Error: No diversity metrics calculated")
        return
    
    # Part 2: Integrate with metadata
    print("\n2. Integrating with metadata...")
    combined_df = integrate_with_metadata(diversity_df)
    
    # Part 3: Statistical analysis
    print("\n3. Performing statistical analysis...")
    stats_results = perform_group_comparisons(combined_df)
    
    # Part 4: Gini coefficient analysis
    print("\n4. Calculating Gini coefficients...")
    gini_df = calculate_all_gini_coefficients()
    
    if gini_df is None:
        print("Error: No Gini coefficients calculated")
        return
    
    # Part 5: Create statistical visualizations (with annotations)
    print("\n5. Creating statistical visualizations with annotations...")
    
    # Create diversity plots with statistical annotations
    print("  - Creating diversity plots with statistical annotations...")
    create_diversity_plots_with_stats_pdf(combined_df)
    
    # Create enhanced clonality plot
    print("  - Creating enhanced clonality plot...")
    create_enhanced_clonality_plot_pdf(combined_df)
    
    # Create individual chain plots
    print("  - Creating individual chain plots...")
    create_individual_chain_plots_pdf_comprehensive(combined_df)
    
    # Create statistical summary CSV
    print("  - Creating statistical significance summary...")
    create_statistical_summary_csv(combined_df)
    
    # Part 6: Create all other visualizations in PDF format
    print("\n6. Creating additional visualizations in PDF format...")
    
    # Diversity plots
    print("  - Creating basic diversity plots...")
    create_diversity_plots_pdf(combined_df)
    
    # Gini coefficient plots
    print("  - Creating Gini coefficient plots...")
    create_gini_plots_pdf(gini_df)
    print("  - Creating B-cell/T-cell comparison...")
    create_bcell_tcell_comparison_pdf(gini_df)
    
    # Additional plots
    print("  - Creating comprehensive clonality plot...")
    create_enhanced_clonality_pdf(combined_df)
    print("  - Creating individual chain plots (Gini)...")
    create_individual_chain_plots_pdf(gini_df)
    print("  - Creating statistical summary...")
    create_statistical_summary_pdf(gini_df)
    
    # Part 7: Additional statistical analysis and save results
    print("\n7. Performing additional statistical analysis...")
    
    # Perform pairwise Mann-Whitney tests for Gini coefficients
    logger = StatisticalResultsLogger()
    
    # Add Gini coefficient summary statistics
    logger.add_section("GINI COEFFICIENT SUMMARY STATISTICS")
    summary_stats = gini_df.groupby(['group', 'chain'])['expansion_gini'].agg(['mean', 'std', 'count']).round(3)
    logger.add_table("Summary Statistics by Group and Chain", summary_stats)
    
    # Perform pairwise comparisons
    logger.add_section("PAIRWISE MANN-WHITNEY TESTS FOR GINI COEFFICIENTS")
    
    all_pairwise_results = []
    chains = gini_df['chain'].unique()
    groups = ['Control', 'Classical_Dengue_Fever', 
              'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        chain_results = []
        
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                group1 = groups[i]
                group2 = groups[j]
                
                data1 = chain_data[chain_data['group'] == group1]['expansion_gini'].dropna()
                data2 = chain_data[chain_data['group'] == group2]['expansion_gini'].dropna()
                
                if len(data1) > 1 and len(data2) > 1:
                    u_stat, p_value = mannwhitneyu(data1, data2)
                    
                    # Calculate effect size
                    mean1, mean2 = data1.mean(), data2.mean()
                    std1, std2 = data1.std(), data2.std()
                    n1, n2 = len(data1), len(data2)
                    pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
                    
                    if pooled_std > 0:
                        cohens_d = (mean2 - mean1) / pooled_std
                    else:
                        cohens_d = np.nan
                    
                    # Determine significance
                    if p_value < 0.001:
                        significance = '***'
                    elif p_value < 0.01:
                        significance = '**'
                    elif p_value < 0.05:
                        significance = '*'
                    else:
                        significance = 'ns'
                    
                    result = {
                        'chain': chain,
                        'comparison': f'{group1}_vs_{group2}',
                        'u_statistic': u_stat,
                        'p_value': p_value,
                        'significance': significance,
                        'cohens_d': cohens_d,
                        'mean_diff': mean2 - mean1
                    }
                    
                    chain_results.append(result)
                    all_pairwise_results.append(result)
        
        if chain_results:
            chain_df = pd.DataFrame(chain_results)
            logger.add_table(f"{chain} Chain Pairwise Comparisons", chain_df)
    
    # Save pairwise results to CSV
    if all_pairwise_results:
        pairwise_df = pd.DataFrame(all_pairwise_results)
        output_path = os.path.join(RESULTS_DIR, "pairwise_mann_whitney_results.csv")
        pairwise_df.to_csv(output_path, index=False)
        print(f"Pairwise Mann-Whitney results saved to {output_path}")
    
    # Calculate effect sizes
    logger.add_section("EFFECT SIZE ANALYSIS")
    effect_sizes = []
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        control_data = chain_data[chain_data['group'] == 'Control']['expansion_gini'].dropna()
        
        for disease_group in ['Classical_Dengue_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']:
            disease_data = chain_data[chain_data['group'] == disease_group]['expansion_gini'].dropna()
            
            if len(control_data) > 1 and len(disease_data) > 1:
                mean_ctrl = control_data.mean()
                mean_disease = disease_data.mean()
                std_ctrl = control_data.std()
                std_disease = disease_data.std()
                n_ctrl = len(control_data)
                n_disease = len(disease_data)
                
                pooled_std = np.sqrt(((n_ctrl-1)*std_ctrl**2 + (n_disease-1)*std_disease**2) / (n_ctrl+n_disease-2))
                
                if pooled_std > 0:
                    cohens_d = (mean_disease - mean_ctrl) / pooled_std
                    
                    effect_sizes.append({
                        'chain': chain,
                        'comparison': f'Control_vs_{disease_group}',
                        "cohens_d": cohens_d,
                        'mean_control': mean_ctrl,
                        f'mean_{disease_group}': mean_disease,
                        'difference': mean_disease - mean_ctrl
                    })
    
    if effect_sizes:
        effect_df = pd.DataFrame(effect_sizes)
        logger.add_table("Effect Size Analysis (Cohen's d)", effect_df)
        
        output_path = os.path.join(RESULTS_DIR, "effect_size_analysis.csv")
        effect_df.to_csv(output_path, index=False)
        print(f"Effect size analysis saved to {output_path}")
    
    # Save the statistical report
    logger.save_to_file()
    
    # Part 8: Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*70)
    print("\nAll results saved in the 'RESULTS' folder:")
    
    # List all generated files
    result_files = os.listdir(RESULTS_DIR)
    
    print("\nPDF Files:")
    pdf_files = [f for f in result_files if f.endswith('.pdf')]
    for pdf in sorted(pdf_files):
        print(f"  - {pdf}")
    
    print("\nCSV Files:")
    csv_files = [f for f in result_files if f.endswith('.csv')]
    for csv in sorted(csv_files):
        print(f"  - {csv}")
    
    print("\nText Reports:")
    txt_files = [f for f in result_files if f.endswith('.txt')]
    for txt in sorted(txt_files):
        print(f"  - {txt}")
    
    print(f"\nTotal files generated: {len(result_files)}")
    print(f"PDF files: {len(pdf_files)}")
    print(f"CSV files: {len(csv_files)}")
    print(f"Text files: {len(txt_files)}")
    
    print("\n" + "="*70)
    print("KEY OUTPUT FILES:")
    print("="*70)
    print("1. dengue_repertoire_diversity_with_stats.pdf - Main diversity plots with statistical annotations")
    print("2. dengue_repertoire_diversity.pdf - Basic diversity plots")
    print("3. dengue_gini_coefficient_analysis.pdf - Gini coefficient analysis")
    print("4. statistical_significance_summary.csv - Summary of all statistical tests")
    print("5. gini_statistical_analysis_report.txt - Detailed statistical report")

# Run the main analysis
if __name__ == "__main__":
    main()
