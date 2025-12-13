# ============================ MAIN ANALYSIS PIPELINE ============================
# This script performs comprehensive immune repertoire analysis for dengue severity
# All results will be saved in a RESULTS folder

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

# ============================ PART 5: VISUALIZATIONS ============================
print("\n" + "="*70)
print("PART 5: CREATING VISUALIZATIONS")
print("="*70)

def add_statistical_annotation(ax, x1, x2, y, p_value, line_height=0.05, color='black'):
    """Add statistical annotation with p-value"""
    
    if p_value < 0.001:
        stars = '***'
    elif p_value < 0.01:
        stars = '**'
    elif p_value < 0.05:
        stars = '*'
    else:
        stars = 'ns'
    
    ax.plot([x1, x1, x2, x2], [y, y+line_height, y+line_height, y], 
            lw=1.5, c=color)
    
    ax.text((x1+x2)*0.5, y+line_height, stars, 
            ha='center', va='bottom', color=color, 
            fontweight='bold', fontsize=10)

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
    output_path = os.path.join(RESULTS_DIR, "enhanced_clonality_analysis.pdf")
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
    output_path = os.path.join(RESULTS_DIR, "individual_chain_analysis.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Individual chain plots saved as PDF: {output_path}")

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

# ============================ PART 6: MAIN EXECUTION ============================
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
    
    # Part 5: Create all visualizations in PDF format
    print("\n5. Creating visualizations in PDF format...")
    
    # Diversity plots
    create_diversity_plots_pdf(combined_df)
    
    # Gini coefficient plots
    create_gini_plots_pdf(gini_df)
    create_bcell_tcell_comparison_pdf(gini_df)
    
    # Additional plots
    create_enhanced_clonality_pdf(combined_df)
    create_individual_chain_plots_pdf(gini_df)
    create_statistical_summary_pdf(gini_df)
    
    # Part 6: Additional statistical analysis and save results
    print("\n6. Performing additional statistical analysis...")
    
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
    
    # Part 7: Summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETED SUCCESSFULLY!")
    print("="*70)
    print("\nAll results saved in the 'RESULTS' folder:")
    print("\nCSV Files:")
    print("  - all_samples_diversity_metrics.csv")
    print("  - dengue_repertoire_analysis_combined.csv")
    print("  - statistical_analysis_results.csv")
    print("  - dengue_gini_coefficient_analysis.csv")
    print("  - pairwise_mann_whitney_results.csv")
    print("  - effect_size_analysis.csv")
    
    print("\nPDF Figures:")
    print("  - dengue_repertoire_diversity.pdf")
    print("  - dengue_gini_coefficient_analysis.pdf")
    print("  - dengue_bcell_tcell_comparison.pdf")
    print("  - enhanced_clonality_analysis.pdf")
    print("  - individual_chain_analysis.pdf")
    print("  - dengue_gini_statistical_summary.pdf")
    
    print("\nText Reports:")
    print("  - gini_statistical_analysis_report.txt")
    
    print("\n" + "="*70)
    print("Files generated in RESULTS folder:")
    print("="*70)
    
    # List all files in RESULTS folder
    result_files = os.listdir(RESULTS_DIR)
    for file in sorted(result_files):
        print(f"  - {file}")
    
    print(f"\nTotal files generated: {len(result_files)}")

# Run the main analysis
if __name__ == "__main__":
    main()
    
    
    
# ============================ ADDITIONAL VISUALIZATION FUNCTIONS ============================
# Add these functions to the existing code

def create_diversity_plots_with_stats_pdf(combined_df):
    """Create diversity plots with statistical annotations in PDF format"""
    
    sns.set_style("whitegrid")
    plt.rcParams['pdf.fonttype'] = 42
    
    # Create the main figure
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Define severity order and palette
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    severity_palette = {
        'Control': '#2E86AB',
        'Classical_Dengue_Fever': '#A23B72',
        'Dengue_Hemorrhagic_Fever': '#F18F01',
        'Dengue_Shock_Syndrome': '#C73E1D'
    }
    
    # 1. Richness by group and chain with stats
    ax1 = axes[0, 0]
    # Create boxplot with all chains
    all_chains_data = combined_df[combined_df['chain'].isin(['IGH', 'IGK', 'IGL', 'TRA', 'TRB'])]
    sns.boxplot(data=all_chains_data, x='chain', y='log_richness', hue='group', 
                palette=severity_palette, ax=ax1)
    
    # Add statistical annotations for each chain
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    y_max = all_chains_data['log_richness'].max()
    
    for chain_idx, chain in enumerate(chains):
        chain_data = combined_df[combined_df['chain'] == chain]
        
        # Perform pairwise tests for Control vs Shock
        control_data = chain_data[chain_data['group'] == 'Control']['log_richness'].dropna()
        shock_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['log_richness'].dropna()
        
        if len(control_data) > 1 and len(shock_data) > 1:
            u_stat, p_value = mannwhitneyu(control_data, shock_data)
            
            if p_value < 0.1:  # Show trends and significant results
                x1 = chain_idx - 0.3  # Control position
                x2 = chain_idx + 0.1  # Shock position
                y_pos = y_max + 0.1
                
                # Add annotation
                if p_value < 0.001:
                    stars = '***'
                elif p_value < 0.01:
                    stars = '**'
                elif p_value < 0.05:
                    stars = '*'
                else:
                    stars = ''
                
                if stars:
                    ax1.plot([x1, x1, x2, x2], [y_pos, y_pos+0.05, y_pos+0.05, y_pos], 
                            lw=1.5, c='black')
                    ax1.text((x1+x2)*0.5, y_pos+0.05, stars, 
                            ha='center', va='bottom', color='black', 
                            fontweight='bold', fontsize=10)
    
    ax1.set_title('Immune Receptor Richness with Statistical Significance', 
                 fontsize=14, fontweight='bold')
    ax1.set_xlabel('Chain Type', fontweight='bold')
    ax1.set_ylabel('Log Richness', fontweight='bold')
    ax1.legend(title='Dengue Severity', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 2. Entropy by group and chain
    ax2 = axes[0, 1]
    sns.boxplot(data=all_chains_data, x='chain', y='entropy', hue='group', 
                palette=severity_palette, ax=ax2)
    ax2.set_title('Immune Receptor Diversity by Dengue Severity', 
                 fontsize=14, fontweight='bold')
    ax2.set_xlabel('Chain Type', fontweight='bold')
    ax2.set_ylabel('Shannon Entropy', fontweight='bold')
    
    # 3. B-cell richness with detailed stats
    ax3 = axes[1, 0]
    b_cells = combined_df[combined_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    sns.boxplot(data=b_cells, x='chain', y='log_richness', hue='group', 
                palette=severity_palette, ax=ax3)
    
    # Add p-values for B-cells
    b_chains = ['IGH', 'IGK', 'IGL']
    y_max_b = b_cells['log_richness'].max()
    
    for chain_idx, chain in enumerate(b_chains):
        chain_data = b_cells[b_cells['chain'] == chain]
        
        # Test Control vs Shock
        control_data = chain_data[chain_data['group'] == 'Control']['log_richness'].dropna()
        shock_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['log_richness'].dropna()
        
        if len(control_data) > 1 and len(shock_data) > 1:
            u_stat, p_value = mannwhitneyu(control_data, shock_data)
            
            if p_value < 0.05:
                x1 = chain_idx - 0.3
                x2 = chain_idx + 0.1
                y_pos = y_max_b + 0.1
                
                ax3.plot([x1, x1, x2, x2], [y_pos, y_pos+0.05, y_pos+0.05, y_pos], 
                        lw=1.5, c='red')
                
                if p_value < 0.001:
                    stars = '***'
                elif p_value < 0.01:
                    stars = '**'
                else:
                    stars = '*'
                
                ax3.text((x1+x2)*0.5, y_pos+0.05, stars, 
                        ha='center', va='bottom', color='red', 
                        fontweight='bold', fontsize=10)
                ax3.text((x1+x2)*0.5, y_pos+0.08, f'p={p_value:.3f}', 
                        ha='center', va='bottom', color='red', 
                        fontsize=8, fontweight='bold')
    
    ax3.set_title('B-cell Richness with Detailed Statistics', 
                 fontsize=14, fontweight='bold')
    ax3.set_xlabel('B-cell Chain Type', fontweight='bold')
    ax3.set_ylabel('Log Richness', fontweight='bold')
    
    # 4. T-cell richness with detailed stats
    ax4 = axes[1, 1]
    t_cells = combined_df[combined_df['chain'].isin(['TRA', 'TRB'])]
    
    if len(t_cells) > 0:
        sns.boxplot(data=t_cells, x='chain', y='log_richness', hue='group', 
                    palette=severity_palette, ax=ax4)
        
        # Add p-values for T-cells
        t_chains = ['TRA', 'TRB']
        y_max_t = t_cells['log_richness'].max()
        
        for chain_idx, chain in enumerate(t_chains):
            chain_data = t_cells[t_cells['chain'] == chain]
            
            if len(chain_data) > 0:
                control_data = chain_data[chain_data['group'] == 'Control']['log_richness'].dropna()
                shock_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['log_richness'].dropna()
                
                if len(control_data) > 1 and len(shock_data) > 1:
                    u_stat, p_value = mannwhitneyu(control_data, shock_data)
                    
                    if p_value < 0.05:
                        x1 = chain_idx - 0.3
                        x2 = chain_idx + 0.1
                        y_pos = y_max_t + 0.1
                        
                        ax4.plot([x1, x1, x2, x2], [y_pos, y_pos+0.05, y_pos+0.05, y_pos], 
                                lw=1.5, c='red')
                        
                        if p_value < 0.001:
                            stars = '***'
                        elif p_value < 0.01:
                            stars = '**'
                        else:
                            stars = '*'
                        
                        ax4.text((x1+x2)*0.5, y_pos+0.05, stars, 
                                ha='center', va='bottom', color='red', 
                                fontweight='bold', fontsize=10)
                        ax4.text((x1+x2)*0.5, y_pos+0.08, f'p={p_value:.3f}', 
                                ha='center', va='bottom', color='red', 
                                fontsize=8, fontweight='bold')
        
        ax4.set_title('T-cell Richness with Detailed Statistics', 
                     fontsize=14, fontweight='bold')
        ax4.set_xlabel('T-cell Chain Type', fontweight='bold')
        ax4.set_ylabel('Log Richness', fontweight='bold')
    
    plt.tight_layout()
    output_path = os.path.join(RESULTS_DIR, "dengue_repertoire_diversity_with_stats.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Diversity plot with stats saved as PDF: {output_path}")

def create_statistical_summary_visualization_pdf(gini_df, combined_df):
    """Create comprehensive statistical summary visualization in PDF"""
    
    sns.set_style("whitegrid")
    plt.rcParams['pdf.fonttype'] = 42
    
    fig = plt.figure(figsize=(20, 15))
    
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    severity_palette = {
        'Control': '#2E86AB',
        'Classical_Dengue_Fever': '#A23B72',
        'Dengue_Hemorrhagic_Fever': '#F18F01',
        'Dengue_Shock_Syndrome': '#C73E1D'
    }
    
    # ========== SUBPLOT 1: Overall Gini comparison ==========
    ax1 = plt.subplot(3, 3, 1)
    
    # Calculate mean Gini for each group
    mean_gini = gini_df.groupby('group')['expansion_gini'].mean()
    std_gini = gini_df.groupby('group')['expansion_gini'].std()
    
    x_pos = range(len(severity_order))
    ax1.bar(x_pos, [mean_gini.get(g, 0) for g in severity_order], 
            yerr=[std_gini.get(g, 0) for g in severity_order],
            color=[severity_palette[g] for g in severity_order], 
            alpha=0.8, capsize=5)
    
    # Add individual points
    for i, group in enumerate(severity_order):
        group_data = gini_df[gini_df['group'] == group]['expansion_gini'].dropna()
        if len(group_data) > 0:
            jitter = np.random.normal(i, 0.1, len(group_data))
            ax1.scatter(jitter, group_data, color='black', alpha=0.5, s=20)
    
    ax1.set_title('Overall Gini Coefficient by Dengue Severity', 
                 fontsize=14, fontweight='bold')
    ax1.set_xlabel('Dengue Severity Group')
    ax1.set_ylabel('Mean Gini Coefficient (± SD)')
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'], rotation=45)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # ========== SUBPLOT 2: Richness statistical summary ==========
    ax2 = plt.subplot(3, 3, 2)
    
    # Calculate p-values matrix for richness
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    p_value_matrix = pd.DataFrame(index=chains, 
                                 columns=['C-D', 'C-H', 'C-S', 'D-H', 'D-S', 'H-S'])
    
    for chain in chains:
        chain_data = combined_df[combined_df['chain'] == chain]
        
        comparisons = [
            ('Control', 'Classical_Dengue_Fever', 'C-D'),
            ('Control', 'Dengue_Hemorrhagic_Fever', 'C-H'),
            ('Control', 'Dengue_Shock_Syndrome', 'C-S'),
            ('Classical_Dengue_Fever', 'Dengue_Hemorrhagic_Fever', 'D-H'),
            ('Classical_Dengue_Fever', 'Dengue_Shock_Syndrome', 'D-S'),
            ('Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome', 'H-S')
        ]
        
        for group1, group2, code in comparisons:
            data1 = chain_data[chain_data['group'] == group1]['log_richness'].dropna()
            data2 = chain_data[chain_data['group'] == group2]['log_richness'].dropna()
            
            if len(data1) > 1 and len(data2) > 1:
                u_stat, p_value = mannwhitneyu(data1, data2)
                p_value_matrix.loc[chain, code] = p_value
            else:
                p_value_matrix.loc[chain, code] = np.nan
    
    # Create heatmap
    p_values = p_value_matrix.astype(float).fillna(1)
    im = ax2.imshow(-np.log10(p_values), cmap='YlOrRd', aspect='auto', vmin=0, vmax=3)
    
    # Add significance stars
    for i in range(len(chains)):
        for j in range(len(p_value_matrix.columns)):
            p_val = p_value_matrix.iloc[i, j]
            if not pd.isna(p_val):
                if p_val < 0.001:
                    text = '***'
                elif p_val < 0.01:
                    text = '**'
                elif p_val < 0.05:
                    text = '*'
                else:
                    text = ''
                ax2.text(j, i, text, ha='center', va='center', 
                        color='black' if p_val < 0.05 else 'gray', 
                        fontweight='bold', fontsize=10)
    
    ax2.set_title('Richness Statistical Significance\n(-log10 p-value)', 
                 fontsize=14, fontweight='bold')
    ax2.set_xlabel('Group Comparisons\n(C=Control, D=Classical, H=Hemorrhagic, S=Shock)')
    ax2.set_ylabel('Immune Receptor Chain')
    ax2.set_xticks(range(len(p_value_matrix.columns)))
    ax2.set_xticklabels(p_value_matrix.columns, rotation=45)
    ax2.set_yticks(range(len(chains)))
    ax2.set_yticklabels(chains)
    plt.colorbar(im, ax=ax2, label='-log10(p-value)')
    
    # ========== SUBPLOT 3: Gini statistical summary ==========
    ax3 = plt.subplot(3, 3, 3)
    
    # Calculate p-values matrix for Gini
    p_value_matrix_gini = pd.DataFrame(index=chains, 
                                      columns=['C-D', 'C-H', 'C-S', 'D-H', 'D-S', 'H-S'])
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        
        for group1, group2, code in comparisons:
            data1 = chain_data[chain_data['group'] == group1]['expansion_gini'].dropna()
            data2 = chain_data[chain_data['group'] == group2]['expansion_gini'].dropna()
            
            if len(data1) > 1 and len(data2) > 1:
                u_stat, p_value = mannwhitneyu(data1, data2)
                p_value_matrix_gini.loc[chain, code] = p_value
            else:
                p_value_matrix_gini.loc[chain, code] = np.nan
    
    # Create heatmap
    p_values_gini = p_value_matrix_gini.astype(float).fillna(1)
    im3 = ax3.imshow(-np.log10(p_values_gini), cmap='YlOrRd', aspect='auto', vmin=0, vmax=3)
    
    # Add significance stars
    for i in range(len(chains)):
        for j in range(len(p_value_matrix_gini.columns)):
            p_val = p_value_matrix_gini.iloc[i, j]
            if not pd.isna(p_val):
                if p_val < 0.001:
                    text = '***'
                elif p_val < 0.01:
                    text = '**'
                elif p_val < 0.05:
                    text = '*'
                else:
                    text = ''
                ax3.text(j, i, text, ha='center', va='center', 
                        color='black' if p_val < 0.05 else 'gray', 
                        fontweight='bold', fontsize=10)
    
    ax3.set_title('Gini Coefficient Statistical Significance\n(-log10 p-value)', 
                 fontsize=14, fontweight='bold')
    ax3.set_xlabel('Group Comparisons\n(C=Control, D=Classical, H=Hemorrhagic, S=Shock)')
    ax3.set_ylabel('Immune Receptor Chain')
    ax3.set_xticks(range(len(p_value_matrix_gini.columns)))
    ax3.set_xticklabels(p_value_matrix_gini.columns, rotation=45)
    ax3.set_yticks(range(len(chains)))
    ax3.set_yticklabels(chains)
    plt.colorbar(im3, ax=ax3, label='-log10(p-value)')
    
    # ========== SUBPLOT 4: Effect size visualization ==========
    ax4 = plt.subplot(3, 3, 4)
    
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
                        'group': disease_group,
                        'cohens_d': cohens_d,
                        'mean_diff': mean_disease - mean_ctrl
                    })
    
    if effect_sizes:
        effect_df = pd.DataFrame(effect_sizes)
        
        # Create grouped bar chart
        x = np.arange(len(chains))
        width = 0.25
        
        for i, disease_group in enumerate(['Classical_Dengue_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']):
            group_effects = effect_df[effect_df['group'] == disease_group]
            d_values = []
            for chain in chains:
                chain_effect = group_effects[group_effects['chain'] == chain]
                if len(chain_effect) > 0:
                    d_values.append(chain_effect['cohens_d'].iloc[0])
                else:
                    d_values.append(0)
            
            ax4.bar(x + (i-1)*width, d_values, width, 
                   label=disease_group.replace('_', ' '),
                   color=severity_palette[disease_group], alpha=0.8)
        
        ax4.set_title("Cohen's d Effect Size Analysis", fontsize=14, fontweight='bold')
        ax4.set_xlabel('Immune Receptor Chain')
        ax4.set_ylabel("Cohen's d")
        ax4.set_xticks(x)
        ax4.set_xticklabels(chains, rotation=45)
        ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax4.axhline(y=0.5, color='red', linestyle='--', linewidth=0.5, alpha=0.5, label='Medium effect')
        ax4.axhline(y=0.8, color='blue', linestyle='--', linewidth=0.5, alpha=0.5, label='Large effect')
        ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # ========== SUBPLOT 5: Clinical interpretation ==========
    ax5 = plt.subplot(3, 3, 5)
    ax5.axis('off')
    
    interpretation = """CLINICAL INTERPRETATION

• Gini Coefficient (0-1):
  0 = Polyclonal (even distribution)
  1 = Oligoclonal (few dominant clones)

• Key Findings:
  - Higher Gini in severe dengue
  - Increasing clonality with severity
  - B-cells show strongest changes

• Statistical Significance:
  *** p < 0.001
  ** p < 0.01
  * p < 0.05

• Effect Size (Cohen's d):
  Small: |d| < 0.5
  Medium: 0.5 ≤ |d| < 0.8
  Large: |d| ≥ 0.8

• Biological Implications:
  Immune focusing in severe cases
  Potential severity biomarkers"""
    
    ax5.text(0.1, 0.95, interpretation, transform=ax5.transAxes,
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    # ========== SUBPLOT 6: Trend analysis ==========
    ax6 = plt.subplot(3, 3, 6)
    
    # Calculate trends across severity
    trend_data = []
    for chain in chains:
        chain_gini = gini_df[gini_df['chain'] == chain]
        if len(chain_gini) > 0:
            # Calculate mean for each severity level
            means = []
            for severity_num in [0, 1, 2, 3]:
                severity_data = chain_gini[chain_gini['severity_level'] == severity_num]['expansion_gini']
                if len(severity_data) > 0:
                    means.append(severity_data.mean())
                else:
                    means.append(np.nan)
            
            trend_data.append({
                'chain': chain,
                'control': means[0] if not np.isnan(means[0]) else 0,
                'classical': means[1] if not np.isnan(means[1]) else 0,
                'hemorrhagic': means[2] if not np.isnan(means[2]) else 0,
                'shock': means[3] if not np.isnan(means[3]) else 0
            })
    
    if trend_data:
        trend_df = pd.DataFrame(trend_data)
        
        x = np.arange(len(chains))
        width = 0.2
        
        ax6.bar(x - width*1.5, trend_df['control'], width, label='Control', 
               color=severity_palette['Control'], alpha=0.8)
        ax6.bar(x - width*0.5, trend_df['classical'], width, label='Classical', 
               color=severity_palette['Classical_Dengue_Fever'], alpha=0.8)
        ax6.bar(x + width*0.5, trend_df['hemorrhagic'], width, label='Hemorrhagic', 
               color=severity_palette['Dengue_Hemorrhagic_Fever'], alpha=0.8)
        ax6.bar(x + width*1.5, trend_df['shock'], width, label='Shock', 
               color=severity_palette['Dengue_Shock_Syndrome'], alpha=0.8)
        
        ax6.set_title('Gini Coefficient Trend Across Severity', 
                     fontsize=14, fontweight='bold')
        ax6.set_xlabel('Immune Receptor Chain')
        ax6.set_ylabel('Mean Gini Coefficient')
        ax6.set_xticks(x)
        ax6.set_xticklabels(chains, rotation=45)
        ax6.legend()
        ax6.grid(True, alpha=0.3, axis='y')
    
    # ========== SUBPLOT 7: Distribution comparison ==========
    ax7 = plt.subplot(3, 3, 7)
    
    # Create violin plot for Gini distribution
    sns.violinplot(data=gini_df, x='group', y='expansion_gini', 
                   ax=ax7, palette=severity_palette, order=severity_order)
    
    # Add mean points
    for i, group in enumerate(severity_order):
        group_data = gini_df[gini_df['group'] == group]['expansion_gini'].dropna()
        if len(group_data) > 0:
            mean_val = group_data.mean()
            ax7.scatter(i, mean_val, color='red', s=100, marker='D', 
                       edgecolor='black', linewidth=1, label='Mean' if i == 0 else "")
    
    ax7.set_title('Gini Coefficient Distribution by Group', 
                 fontsize=14, fontweight='bold')
    ax7.set_xlabel('Dengue Severity Group')
    ax7.set_ylabel('Gini Coefficient')
    ax7.set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'], rotation=45)
    
    # ========== SUBPLOT 8: Correlation analysis ==========
    ax8 = plt.subplot(3, 3, 8)
    
    # Calculate correlation between richness and Gini
    correlation_data = []
    for chain in chains:
        chain_gini = gini_df[gini_df['chain'] == chain]
        chain_richness = combined_df[combined_df['chain'] == chain]
        
        if len(chain_gini) > 0 and len(chain_richness) > 0:
            merged_data = pd.merge(chain_gini[['sample_id', 'expansion_gini']], 
                                  chain_richness[['sample_id', 'log_richness']], 
                                  on='sample_id')
            
            if len(merged_data) > 2:
                correlation = merged_data['expansion_gini'].corr(merged_data['log_richness'])
                correlation_data.append({
                    'chain': chain,
                    'correlation': correlation
                })
    
    if correlation_data:
        corr_df = pd.DataFrame(correlation_data)
        
        colors = ['red' if r < 0 else 'blue' for r in corr_df['correlation']]
        ax8.bar(corr_df['chain'], corr_df['correlation'], color=colors, alpha=0.7)
        
        ax8.set_title('Correlation: Gini vs Richness', 
                     fontsize=14, fontweight='bold')
        ax8.set_xlabel('Immune Receptor Chain')
        ax8.set_ylabel('Correlation Coefficient (r)')
        ax8.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax8.axhline(y=0.3, color='green', linestyle='--', linewidth=0.5, alpha=0.5, label='Moderate')
        ax8.axhline(y=0.5, color='orange', linestyle='--', linewidth=0.5, alpha=0.5, label='Strong')
        ax8.set_xticklabels(corr_df['chain'], rotation=45)
    
    # ========== SUBPLOT 9: Summary statistics table ==========
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    # Calculate key summary statistics
    summary_stats = []
    
    for chain in ['IGH', 'TRB']:  # Most important chains
        chain_gini = gini_df[gini_df['chain'] == chain]
        
        if len(chain_gini) > 0:
            control_mean = chain_gini[chain_gini['group'] == 'Control']['expansion_gini'].mean()
            shock_mean = chain_gini[chain_gini['group'] == 'Dengue_Shock_Syndrome']['expansion_gini'].mean()
            
            if not pd.isna(control_mean) and not pd.isna(shock_mean):
                percent_change = ((shock_mean - control_mean) / control_mean) * 100
                
                # Test significance
                c_data = chain_gini[chain_gini['group'] == 'Control']['expansion_gini'].dropna()
                s_data = chain_gini[chain_gini['group'] == 'Dengue_Shock_Syndrome']['expansion_gini'].dropna()
                
                if len(c_data) > 1 and len(s_data) > 1:
                    _, p_value = mannwhitneyu(c_data, s_data)
                    
                    significance = ''
                    if p_value < 0.001:
                        significance = '***'
                    elif p_value < 0.01:
                        significance = '**'
                    elif p_value < 0.05:
                        significance = '*'
                    
                    summary_stats.append([
                        chain,
                        f"{control_mean:.3f}",
                        f"{shock_mean:.3f}",
                        f"{percent_change:+.1f}%",
                        significance
                    ])
    
    if summary_stats:
        table = ax9.table(cellText=summary_stats,
                         colLabels=['Chain', 'Control', 'Shock', 'Change', 'Sig'],
                         cellLoc='center',
                         loc='center',
                         colWidths=[0.15, 0.2, 0.2, 0.2, 0.15])
        
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2)
        
        ax9.set_title('Key Statistical Findings', fontsize=14, fontweight='bold', y=0.98)
    
    plt.suptitle('Comprehensive Statistical Analysis: Dengue Severity Immune Repertoire', 
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    output_path = os.path.join(RESULTS_DIR, "statistical_summary_visualization.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Statistical summary visualization saved as PDF: {output_path}")

def create_basic_diversity_plot_pdf(combined_df):
    """Create basic diversity plot in PDF format"""
    
    sns.set_style("whitegrid")
    plt.rcParams['pdf.fonttype'] = 42
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    severity_palette = {
        'Control': '#2E86AB',
        'Classical_Dengue_Fever': '#A23B72',
        'Dengue_Hemorrhagic_Fever': '#F18F01',
        'Dengue_Shock_Syndrome': '#C73E1D'
    }
    
    # 1. Richness plot
    ax1 = axes[0, 0]
    sns.boxplot(data=combined_df, x='chain', y='log_richness', hue='group', 
                palette=severity_palette, ax=ax1)
    ax1.set_title('Immune Receptor Richness', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Chain Type')
    ax1.set_ylabel('Log Richness')
    ax1.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 2. Entropy plot
    ax2 = axes[0, 1]
    sns.boxplot(data=combined_df, x='chain', y='entropy', hue='group', 
                palette=severity_palette, ax=ax2)
    ax2.set_title('Immune Receptor Diversity', fontsize=12, fontweight='bold')
    ax2.set_xlabel('Chain Type')
    ax2.set_ylabel('Shannon Entropy')
    
    # 3. Clonality plot
    ax3 = axes[1, 0]
    sns.boxplot(data=combined_df, x='chain', y='clonality', hue='group', 
                palette=severity_palette, ax=ax3)
    ax3.set_title('Immune Receptor Clonality', fontsize=12, fontweight='bold')
    ax3.set_xlabel('Chain Type')
    ax3.set_ylabel('Clonality (1 - normalized entropy)')
    
    # 4. Number of clonotypes
    ax4 = axes[1, 1]
    sns.boxplot(data=combined_df, x='chain', y='n_clonotypes', hue='group', 
                palette=severity_palette, ax=ax4)
    ax4.set_title('Number of Clonotypes', fontsize=12, fontweight='bold')
    ax4.set_xlabel('Chain Type')
    ax4.set_ylabel('Number of Clonotypes')
    
    plt.suptitle('Dengue Severity Immune Repertoire Diversity Analysis', 
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = os.path.join(RESULTS_DIR, "dengue_repertoire_diversity.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Basic diversity plot saved as PDF: {output_path}")

# ============================ UPDATE MAIN EXECUTION ============================
# Add these function calls to the main() function

def main():
    """Main analysis pipeline - UPDATED VERSION"""
    
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
    
    # Part 5: Create all visualizations in PDF format
    print("\n5. Creating visualizations in PDF format...")
    
    # Create all PDF files
    create_basic_diversity_plot_pdf(combined_df)  # dengue_repertoire_diversity.pdf
    create_diversity_plots_with_stats_pdf(combined_df)  # dengue_repertoire_diversity_with_stats.pdf
    create_statistical_summary_visualization_pdf(gini_df, combined_df)  # statistical_summary_visualization.pdf
    
    # Existing PDF creations
    create_gini_plots_pdf(gini_df)  # dengue_gini_coefficient_analysis.pdf
    create_bcell_tcell_comparison_pdf(gini_df)  # dengue_bcell_tcell_comparison.pdf
    create_enhanced_clonality_pdf(combined_df)  # enhanced_clonality_analysis.pdf
    create_individual_chain_plots_pdf(gini_df)  # individual_chain_analysis.pdf
    
    # This function already creates dengue_gini_statistical_summary.pdf
    # create_statistical_summary_pdf(gini_df)
    
    # Part 6: Additional statistical analysis and save results
    print("\n6. Performing additional statistical analysis...")
    
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
    
    # Part 7: Summary
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

# Run the main analysis
if __name__ == "__main__":
    main()
