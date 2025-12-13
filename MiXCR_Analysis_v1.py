import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal, gaussian_kde
import itertools
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Your dengue severity metadata
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

# Create RESULTS directory if it doesn't exist
if not os.path.exists('RESULTS'):
    os.makedirs('RESULTS')
    print("Created RESULTS directory")
    
# Create subdirectories for different types of results
subdirectories = ['raw_data', 'plots', 'statistics', 'tables', 'reports']
for subdir in subdirectories:
    if not os.path.exists(f'RESULTS/{subdir}'):
        os.makedirs(f'RESULTS/{subdir}')
        print(f"Created RESULTS/{subdir} directory")

print("="*60)
print("DENGUE IMMUNE REPERTOIRE ANALYSIS PIPELINE")
print("="*60)
print(f"Analysis started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

class ResultsManager:
    """Manage all results saving with consistent paths"""
    
    @staticmethod
    def save_csv(df, filename, subfolder='tables'):
        """Save DataFrame to CSV in RESULTS folder"""
        path = f'RESULTS/{subfolder}/{filename}'
        df.to_csv(path, index=False)
        print(f"  ✓ Saved: {path}")
        return path
    
    @staticmethod
    def save_plot(fig, filename, subfolder='plots', dpi=300):
        """Save plot to RESULTS folder"""
        path = f'RESULTS/{subfolder}/{filename}'
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f"  ✓ Saved: {path}")
        return path
    
    @staticmethod
    def save_text(content, filename, subfolder='reports'):
        """Save text content to RESULTS folder"""
        path = f'RESULTS/{subfolder}/{filename}'
        with open(path, 'w') as f:
            f.write(content)
        print(f"  ✓ Saved: {path}")
        return path
    
    @staticmethod
    def save_figure(fig, filename, subfolder='plots', dpi=300):
        """Save matplotlib figure"""
        path = f'RESULTS/{subfolder}/{filename}'
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f"  ✓ Saved: {path}")
        return path

# ============================================================================
# PART 1: GINI COEFFICIENT CALCULATION
# ============================================================================

def calculate_gini_coefficient(values):
    """Calculate Gini coefficient for a distribution of values."""
    if len(values) == 0:
        return np.nan
    
    values = np.array(values)
    values = values[values > 0]  # Remove zeros
    n = len(values)
    
    if n == 0:
        return np.nan
    
    # Gini coefficient calculation
    sorted_values = np.sort(values)
    index = np.arange(1, n + 1)
    gini = (np.sum((2 * index - n - 1) * sorted_values)) / (n * np.sum(sorted_values))
    
    return abs(gini)

def extract_gene_name(gene_field):
    """Extract gene name from MiXCR gene field"""
    import re
    if pd.isna(gene_field) or gene_field == "":
        return "None"
    return re.sub(r'\(.*', '', gene_field)

def calculate_clonal_expansion_gini(clone_data):
    """Calculate Gini coefficient for clonal expansion (vertex size)"""
    if len(clone_data) == 0:
        return np.nan
    
    vertex_sizes = clone_data['readCount'].values
    return calculate_gini_coefficient(vertex_sizes)

def process_sample_for_gini(sample_id, chain_file, chain_type):
    """Process a single sample's chain data for Gini coefficient calculation"""
    try:
        df = pd.read_csv(chain_file, sep='\t')
        
        if len(df) == 0:
            return None
        
        processed_data = df[['readCount', 'nSeqCDR3', 
                            'allVHitsWithScore', 'allJHitsWithScore']].copy()
        
        processed_data['v_gene_clean'] = processed_data['allVHitsWithScore'].apply(extract_gene_name)
        processed_data['j_gene_clean'] = processed_data['allJHitsWithScore'].apply(extract_gene_name)
        processed_data['cdr3_nt'] = processed_data['nSeqCDR3'].fillna('')
        
        expansion_gini = calculate_clonal_expansion_gini(processed_data)
        
        return {
            'sample_id': sample_id,
            'chain': chain_type,
            'expansion_gini': expansion_gini,
            'n_clones': len(processed_data),
            'total_reads': processed_data['readCount'].sum()
        }
        
    except Exception as e:
        print(f"Error processing {chain_file}: {e}")
        return None

def calculate_all_gini_coefficients(results_dir="."):
    """Calculate Gini coefficients for all samples and chains"""
    print("\n" + "="*60)
    print("PART 1: CALCULATING GINI COEFFICIENTS")
    print("="*60)
    
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
        
        metadata_df = pd.DataFrame(DENGUE_METADATA)
        results_df = pd.merge(results_df, metadata_df, on='sample_id', how='left')
        
        severity_order = {
            'Control': 0,
            'Classical_Dengue_Fever': 1,
            'Dengue_Hemorrhagic_Fever': 2,
            'Dengue_Shock_Syndrome': 3
        }
        results_df['severity_order'] = results_df['group'].map(severity_order)
        results_df = results_df.sort_values(['severity_order', 'sample_id'])
        
        ResultsManager.save_csv(results_df, 'gini_coefficients.csv', 'raw_data')
        
        print(f"\nGini coefficient calculation completed!")
        print(f"Total records: {len(results_df)}")
        
        return results_df
    else:
        print("No results to save")
        return None

# ============================================================================
# PART 2: STATISTICAL ANALYSIS
# ============================================================================

class StatisticalResultsLogger:
    """Class to save all statistical results to text file"""
    
    def __init__(self):
        self.filename = 'RESULTS/reports/statistical_analysis_report.txt'
        self.results = []
        self.summary_stats = []
        
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = f"""DENGUE SEVERITY GINI COEFFICIENT STATISTICAL ANALYSIS
===========================================================
Analysis Date: {timestamp}
Number of Samples: 24
Severity Groups: Control (n=6), Classical_Dengue_Fever (n=6), 
                 Dengue_Hemorrhagic_Fever (n=6), Dengue_Shock_Syndrome (n=6)
Chains Analyzed: IGH, IGK, IGL, TRA, TRB
Gini Coefficient Range: 0 (polyclonal) to 1 (oligoclonal)
===========================================================

"""
        self.results.append(header)
    
    def add_section(self, title):
        self.results.append(f"\n{'='*60}\n{title}\n{'='*60}\n")
    
    def add_result(self, text):
        self.results.append(text + "\n")
    
    def add_table(self, title, df):
        self.add_section(title)
        self.results.append(df.to_string() + "\n")
    
    def add_summary_stat(self, stat_name, value):
        self.summary_stats.append(f"{stat_name}: {value}")
    
    def save_to_file(self):
        if self.summary_stats:
            self.add_section("ANALYSIS SUMMARY")
            for stat in self.summary_stats:
                self.add_result(stat)
        
        with open(self.filename, 'w') as f:
            f.writelines(self.results)
        
        print(f"  ✓ Saved: {self.filename}")

def perform_statistical_analysis(gini_df):
    """Perform comprehensive statistical analysis"""
    print("\n" + "="*60)
    print("PART 2: PERFORMING STATISTICAL ANALYSIS")
    print("="*60)
    
    logger = StatisticalResultsLogger()
    
    # 1. Descriptive Statistics
    print("1. Calculating descriptive statistics...")
    logger.add_section("1. DESCRIPTIVE STATISTICS")
    
    overall_stats = gini_df['expansion_gini'].describe()
    logger.add_result("Overall Gini Coefficient Statistics:")
    logger.add_result(f"  Mean: {overall_stats['mean']:.4f}")
    logger.add_result(f"  Std: {overall_stats['std']:.4f}")
    logger.add_result(f"  Min: {overall_stats['min']:.4f}")  # FIXED: Changed 'overstalls' to 'overall_stats'
    logger.add_result(f"  25%: {overall_stats['25%']:.4f}")
    logger.add_result(f"  50% (Median): {overall_stats['50%']:.4f}")
    logger.add_result(f"  75%: {overall_stats['75%']:.4f}")
    logger.add_result(f"  Max: {overall_stats['max']:.4f}")
    
    # By group
    logger.add_result("\nBy Dengue Severity Group:")
    group_stats = gini_df.groupby('group')['expansion_gini'].agg(['mean', 'std', 'min', 'max', 'count'])
    for group in ['Control', 'Classical_Dengue_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']:
        if group in group_stats.index:
            stats = group_stats.loc[group]
            logger.add_result(f"  {group}:")
            logger.add_result(f"    Mean: {stats['mean']:.4f} ± {stats['std']:.4f}")
            logger.add_result(f"    Range: {stats['min']:.4f} - {stats['max']:.4f}")
            logger.add_result(f"    N: {int(stats['count'])}")
    
    # By chain
    logger.add_result("\nBy Immune Receptor Chain:")
    chain_stats = gini_df.groupby('chain')['expansion_gini'].agg(['mean', 'std', 'count'])
    for chain in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']:
        if chain in chain_stats.index:
            stats = chain_stats.loc[chain]
            logger.add_result(f"  {chain}: Mean = {stats['mean']:.4f} ± {stats['std']:.4f}, N = {int(stats['count'])}")
    
    # 2. Kruskal-Wallis Tests
    print("2. Performing Kruskal-Wallis tests...")
    logger.add_section("2. KRUSKAL-WALLIS TESTS (Overall Group Differences)")
    logger.add_result("Null Hypothesis: All groups have the same distribution of Gini coefficients")
    logger.add_result("Alternative Hypothesis: At least one group differs")
    logger.add_result("Significance level: α = 0.05\n")
    
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    
    kruskal_results = []
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        
        # Extract data for each group
        group_data = []
        group_labels = []
        for group in severity_order:
            data = chain_data[chain_data['group'] == group]['expansion_gini'].dropna()
            if len(data) > 0:
                group_data.append(data.values)
                group_labels.append(group)
        
        if len(group_data) > 1:
            h_stat, p_value = kruskal(*group_data)
            
            # Determine significance
            if p_value < 0.001:
                significance = "*** (Highly Significant)"
            elif p_value < 0.01:
                significance = "** (Very Significant)"
            elif p_value < 0.05:
                significance = "* (Significant)"
            else:
                significance = "ns (Not Significant)"
            
            kruskal_results.append({
                'Chain': chain,
                'H-Statistic': f"{h_stat:.4f}",
                'p-value': f"{p_value:.6f}",
                'Significance': significance,
                'Interpretation': "Groups differ" if p_value < 0.05 else "No significant difference"
            })
            
            # Log to summary
            if p_value < 0.05:
                logger.add_summary_stat(f"{chain} overall difference", f"p={p_value:.4f} {significance[:3]}")
    
    # Create and display results table
    if kruskal_results:
        results_df = pd.DataFrame(kruskal_results)
        logger.add_table("Kruskal-Wallis Test Results", results_df)
        
        # Save to CSV
        ResultsManager.save_csv(results_df, 'kruskal_wallis_results.csv', 'statistics')
        
        # Count significant findings
        sig_count = sum(1 for r in kruskal_results if r['p-value'] != 'nan' and float(r['p-value']) < 0.05)
        logger.add_result(f"\nSummary: {sig_count} out of {len(kruskal_results)} chains show significant overall differences (p < 0.05)")
    
    # 3. Pairwise Mann-Whitney Tests
    print("3. Performing pairwise Mann-Whitney tests...")
    logger.add_section("3. PAIRWISE MANN-WHITNEY U TESTS")
    logger.add_result("Comparison between each dengue severity group")
    logger.add_result("Non-parametric test for two independent samples")
    logger.add_result("Significance levels: * p<0.05, ** p<0.01, *** p<0.001\n")
    
    pairwise_results = []
    
    for chain in chains:
        logger.add_result(f"\n{chain} Chain:")
        chain_data = gini_df[gini_df['chain'] == chain]
        
        comparisons_made = 0
        significant_comparisons = 0
        
        # Compare each pair of groups
        for i in range(len(severity_order)):
            for j in range(i+1, len(severity_order)):
                group1 = severity_order[i]
                group2 = severity_order[j]
                
                data1 = chain_data[chain_data['group'] == group1]['expansion_gini'].dropna()
                data2 = chain_data[chain_data['group'] == group2]['expansion_gini'].dropna()
                
                if len(data1) > 1 and len(data2) > 1:
                    u_stat, p_value = mannwhitneyu(data1, data2, alternative='two-sided')
                    comparisons_made += 1
                    
                    # Calculate effect size (Cohen's d approximation)
                    mean1, mean2 = data1.mean(), data2.mean()
                    std1, std2 = data1.std(), data2.std()
                    n1, n2 = len(data1), len(data2)
                    pooled_std = np.sqrt(((n1-1)*std1**2 + (n2-1)*std2**2) / (n1+n2-2))
                    
                    if pooled_std > 0:
                        cohens_d = (mean2 - mean1) / pooled_std
                    else:
                        cohens_d = np.nan
                    
                    # Determine significance stars
                    if p_value < 0.001:
                        stars = '***'
                    elif p_value < 0.01:
                        stars = '**'
                    elif p_value < 0.05:
                        stars = '*'
                        significant_comparisons += 1
                    else:
                        stars = 'ns'
                    
                    # Store result
                    result = {
                        'Chain': chain,
                        'Comparison': f"{group1} vs {group2}",
                        'U-Statistic': f"{u_stat:.2f}",
                        'p-value': f"{p_value:.6f}",
                        'Significance': stars,
                        'Effect Size (d)': f"{cohens_d:.3f}",
                        'Mean Diff': f"{mean2 - mean1:.4f}"
                    }
                    
                    pairwise_results.append(result)
                    
                    # Log individual result
                    logger.add_result(f"  {group1} vs {group2}: U={u_stat:.2f}, p={p_value:.4f} {stars}, d={cohens_d:.3f}")
                    
                    # Log to summary if significant
                    if p_value < 0.05:
                        logger.add_summary_stat(f"{chain} {group1}vs{group2}", 
                                              f"p={p_value:.4f}, d={cohens_d:.3f}")
        
        logger.add_result(f"  Total comparisons: {comparisons_made}")
        logger.add_result(f"  Significant comparisons: {significant_comparisons}")
    
    # Save all results to CSV
    if pairwise_results:
        all_results_df = pd.DataFrame(pairwise_results)
        ResultsManager.save_csv(all_results_df, 'pairwise_mann_whitney_results.csv', 'statistics')
        logger.add_result(f"\nDetailed pairwise results saved to 'pairwise_mann_whitney_results.csv'")
    
    # 4. Effect Size Analysis
    print("4. Calculating effect sizes...")
    logger.add_section("4. EFFECT SIZE ANALYSIS")
    logger.add_result("Cohen's d effect size interpretation:")
    logger.add_result("  |d| < 0.2: Negligible effect")
    logger.add_result("  0.2 ≤ |d| < 0.5: Small effect")
    logger.add_result("  0.5 ≤ |d| < 0.8: Medium effect")
    logger.add_result("  |d| ≥ 0.8: Large effect\n")
    
    effect_sizes = []
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        
        # Compare each disease group to control
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
                
                # Pooled standard deviation
                pooled_std = np.sqrt(((n_ctrl-1)*std_ctrl**2 + (n_disease-1)*std_disease**2) / (n_ctrl+n_disease-2))
                
                if pooled_std > 0:
                    cohens_d = (mean_disease - mean_ctrl) / pooled_std
                    
                    # Interpret effect size
                    abs_d = abs(cohens_d)
                    if abs_d < 0.2:
                        magnitude = "Negligible"
                    elif abs_d < 0.5:
                        magnitude = "Small"
                    elif abs_d < 0.8:
                        magnitude = "Medium"
                    else:
                        magnitude = "Large"
                    
                    direction = "increase" if cohens_d > 0 else "decrease"
                    
                    effect_sizes.append({
                        'Chain': chain,
                        'Comparison': f"Control vs {disease_group}",
                        "Cohen's d": f"{cohens_d:.3f}",
                        'Magnitude': magnitude,
                        'Direction': direction,
                        'Mean Control': f"{mean_ctrl:.3f}",
                        f'Mean {disease_group}': f"{mean_disease:.3f}",
                        'Difference': f"{mean_disease - mean_ctrl:.3f}"
                    })
    
    if effect_sizes:
        effect_df = pd.DataFrame(effect_sizes)
        logger.add_table("Effect Size Analysis (Cohen's d)", effect_df)
        
        # Save to CSV
        ResultsManager.save_csv(effect_df, 'effect_size_analysis.csv', 'statistics')
        logger.add_result(f"\nEffect size results saved to 'effect_size_analysis.csv'")
    
    # 5. Trend Analysis
    print("5. Analyzing trends across severity levels...")
    logger.add_section("5. TREND ANALYSIS ACROSS SEVERITY LEVELS")
    logger.add_result("Analysis of Gini coefficient trends from Control → Classical → Hemorrhagic → Shock\n")
    
    trend_results = []
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain].copy()
        severity_mapping = {
            'Control': 0,
            'Classical_Dengue_Fever': 1,
            'Dengue_Hemorrhagic_Fever': 2,
            'Dengue_Shock_Syndrome': 3
        }
        chain_data['severity_numeric'] = chain_data['group'].map(severity_mapping)
        
        if len(chain_data) >= 8:  # Need sufficient data
            # Calculate correlation
            correlation = chain_data['severity_numeric'].corr(chain_data['expansion_gini'])
            
            # Calculate means for each severity level
            means = []
            for severity in range(4):
                level_data = chain_data[chain_data['severity_numeric'] == severity]['expansion_gini']
                if len(level_data) > 0:
                    means.append(level_data.mean())
                else:
                    means.append(np.nan)
            
            # Determine trend
            if not any(np.isnan(means)):
                if all(means[i] <= means[i+1] for i in range(len(means)-1)):
                    trend = "Increasing (↑ with severity)"
                elif all(means[i] >= means[i+1] for i in range(len(means)-1)):
                    trend = "Decreasing (↓ with severity)"
                else:
                    trend = "Non-monotonic"
            else:
                trend = "Incomplete data"
            
            # Interpret correlation
            if abs(correlation) < 0.3:
                strength = "Weak"
            elif abs(correlation) < 0.5:
                strength = "Moderate"
            elif abs(correlation) < 0.7:
                strength = "Strong"
            else:
                strength = "Very strong"
            
            direction = "Positive" if correlation > 0 else "Negative"
            
            trend_results.append({
                'Chain': chain,
                'Correlation (r)': f"{correlation:.3f}",
                'Trend': trend,
                'Strength': f"{strength} {direction}",
                'Control Mean': f"{means[0]:.3f}" if not np.isnan(means[0]) else "NA",
                'Shock Mean': f"{means[3]:.3f}" if not np.isnan(means[3]) else "NA",
                'Relative Change': f"{(means[3]-means[0])/means[0]*100:+.1f}%" if not np.isnan(means[0]) and not np.isnan(means[3]) and means[0] != 0 else "NA"
            })
    
    if trend_results:
        trend_df = pd.DataFrame(trend_results)
        logger.add_table("Trend Analysis Results", trend_df)
        
        # Count increasing trends
        increasing = sum(1 for r in trend_results if "Increasing" in r['Trend'])
        logger.add_result(f"\n{increasing} out of {len(trend_results)} chains show increasing Gini with dengue severity")
    
    # 6. Clinical Interpretation
    print("6. Generating clinical interpretation...")
    logger.add_section("6. CLINICAL INTERPRETATION AND CONCLUSIONS")
    
    # Key findings summary
    logger.add_result("KEY FINDINGS:\n")
    
    # Check for significant Control vs Shock differences
    significant_chains = []
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        c_data = chain_data[chain_data['group'] == 'Control']['expansion_gini'].dropna()
        s_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['expansion_gini'].dropna()
        
        if len(c_data) > 1 and len(s_data) > 1:
            _, p_value = mannwhitneyu(c_data, s_data)
            if p_value < 0.05:
                significant_chains.append(chain)
                mean_ctrl = c_data.mean()
                mean_shock = s_data.mean()
                percent_change = ((mean_shock - mean_ctrl) / mean_ctrl) * 100
                
                logger.add_result(f"• {chain}: Significant difference between Control and Shock Syndrome")
                logger.add_result(f"  Control: {mean_ctrl:.3f}, Shock: {mean_shock:.3f} ({percent_change:+.1f}%)")
                if mean_shock > mean_ctrl:
                    logger.add_result(f"  → Increased clonal expansion in severe dengue")
                else:
                    logger.add_result(f"  → Decreased clonal expansion in severe dengue")
    
    # Biological interpretation
    logger.add_result("\nBIOLOGICAL INTERPRETATION:")
    logger.add_result("• Gini coefficient measures inequality in clone sizes")
    logger.add_result("• Higher Gini (closer to 1) = Oligoclonal expansion")
    logger.add_result("  - Few dominant B/T cell clones")
    logger.add_result("  - Potential antigen-specific immune response")
    logger.add_result("  - Characteristic of focused immune response")
    logger.add_result("• Lower Gini (closer to 0) = Polyclonal repertoire")
    logger.add_result("  - Even distribution of many clones")
    logger.add_result("  - Broad, diverse immune response")
    logger.add_result("  - Characteristic of homeostatic immunity")
    
    # Clinical implications
    logger.add_result("\nCLINICAL IMPLICATIONS FOR DENGUE:")
    if significant_chains:
        logger.add_result(f"• Significant changes in {len(significant_chains)} immune receptor chains")
        logger.add_result("• Severe dengue associated with altered clonal expansion patterns")
        logger.add_result("• Potential biomarkers for disease severity")
        logger.add_result("• May reflect immune system focusing on dengue antigens")
    else:
        logger.add_result("• No significant differences in clonal expansion patterns")
        logger.add_result("• Suggests similar B/T cell repertoire across severity groups")
        logger.add_result("• May indicate non-specific immune activation")
    
    # Recommendations
    logger.add_result("\nRECOMMENDATIONS FOR FURTHER ANALYSIS:")
    logger.add_result("1. Analyze specific expanded clones for dengue antigen specificity")
    logger.add_result("2. Correlate with clinical parameters (viral load, platelet count)")
    logger.add_result("3. Compare with other viral infections for specificity")
    logger.add_result("4. Longitudinal analysis to track repertoire changes")
    logger.add_result("5. Validate findings with larger cohort studies")
    
    # 7. Save the report
    logger.save_to_file()
    
    return logger

# ============================================================================
# PART 3: VISUALIZATION
# ============================================================================

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
            ha='center', va='bottom', color=color, fontweight='bold', fontsize=10)

def create_main_visualizations(gini_df):
    """Create main visualization figures"""
    print("\n" + "="*60)
    print("PART 3: CREATING VISUALIZATIONS")
    print("="*60)
    
    severity_order = ['Control', 'Classical_Dengue_Fever', 
                     'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    severity_palette = {
        'Control': '#2E86AB',
        'Classical_Dengue_Fever': '#A23B72',
        'Dengue_Hemorrhagic_Fever': '#F18F01',
        'Dengue_Shock_Syndrome': '#C73E1D'
    }
    
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    
    # 1. Main statistical summary figure
    print("1. Creating main statistical summary figure...")
    fig1, axes1 = plt.subplots(2, 3, figsize=(15, 10))
    axes1 = axes1.flatten()
    
    for idx, chain in enumerate(chains[:5]):
        ax = axes1[idx]
        chain_data = gini_df[gini_df['chain'] == chain]
        
        if len(chain_data) > 0:
            sns.boxplot(data=chain_data, x='group', y='expansion_gini', hue='group',
                       ax=ax, palette=severity_palette, order=severity_order, legend=False)
            sns.stripplot(data=chain_data, x='group', y='expansion_gini', 
                         ax=ax, color='black', alpha=0.5, size=4, order=severity_order)
            
            # Add statistical annotation for Control vs Shock
            c_data = chain_data[chain_data['group'] == 'Control']['expansion_gini'].dropna()
            s_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['expansion_gini'].dropna()
            
            if len(c_data) > 1 and len(s_data) > 1:
                _, p_value = mannwhitneyu(c_data, s_data)
                if p_value < 0.1:
                    x1 = severity_order.index('Control')
                    x2 = severity_order.index('Dengue_Shock_Syndrome')
                    y_max = chain_data['expansion_gini'].max()
                    y_pos = y_max + 0.1
                    add_statistical_annotation(ax, x1, x2, y_pos, p_value)
            
            ax.set_title(f'{chain} Chain', fontweight='bold')
            ax.set_xlabel('')
            ax.set_ylabel('Gini Coefficient')
            ax.set_xticklabels(['Ctrl', 'Class', 'Hem', 'Shock'], rotation=45)
            ax.grid(True, alpha=0.3, axis='y')
    
    # Summary subplot
    ax_summary = axes1[5]
    ax_summary.axis('off')
    summary_text = """STATISTICAL SUMMARY

Significance Levels:
*** p < 0.001
**  p < 0.01
*   p < 0.05
ns  p ≥ 0.05

Interpretation:
Higher Gini = More oligoclonal
Lower Gini = More polyclonal"""
    
    ax_summary.text(0.1, 0.5, summary_text, transform=ax_summary.transAxes,
                   fontsize=10, verticalalignment='center',
                   bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.suptitle('Gini Coefficient Analysis: Dengue Severity Groups', fontsize=16, fontweight='bold')
    plt.tight_layout()
    ResultsManager.save_figure(fig1, 'statistical_summary.png', 'plots')
    plt.close(fig1)
    
    # 2. Heatmap of p-values
    print("2. Creating p-value heatmap...")
    p_value_data = []
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        row = []
        for i in range(len(severity_order)-1):
            for j in range(i+1, len(severity_order)):
                g1 = severity_order[i]
                g2 = severity_order[j]
                d1 = chain_data[chain_data['group'] == g1]['expansion_gini'].dropna()
                d2 = chain_data[chain_data['group'] == g2]['expansion_gini'].dropna()
                
                if len(d1) > 1 and len(d2) > 1:
                    _, p = mannwhitneyu(d1, d2)
                    row.append(p)
                else:
                    row.append(np.nan)
        p_value_data.append(row)
    
    p_value_df = pd.DataFrame(p_value_data, index=chains, 
                             columns=['C-D', 'C-H', 'C-S', 'D-H', 'D-S', 'H-S'])
    
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    im = ax2.imshow(-np.log10(p_value_df.fillna(1)), cmap='YlOrRd', aspect='auto', vmin=0, vmax=3)
    
    for i in range(len(chains)):
        for j in range(len(p_value_df.columns)):
            p_val = p_value_df.iloc[i, j]
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
                        color='black' if p_val < 0.05 else 'gray', fontweight='bold')
    
    ax2.set_title('Statistical Significance Heatmap\n(-log10 p-values)', fontweight='bold', fontsize=14)
    ax2.set_xlabel('Group Comparisons\n(C=Control, D=Classical, H=Hemorrhagic, S=Shock)', fontweight='bold')
    ax2.set_ylabel('Immune Receptor Chain', fontweight='bold')
    ax2.set_xticks(range(len(p_value_df.columns)))
    ax2.set_xticklabels(p_value_df.columns, rotation=45)
    ax2.set_yticks(range(len(chains)))
    ax2.set_yticklabels(chains)
    plt.colorbar(im, ax=ax2, label='-log10(p-value)')
    plt.tight_layout()
    ResultsManager.save_figure(fig2, 'p_value_heatmap.png', 'plots')
    plt.close(fig2)
    
    # 3. Trend analysis plot
    print("3. Creating trend analysis plot...")
    fig3, axes3 = plt.subplots(1, 2, figsize=(12, 5))
    
    # Scatter plot with trend lines
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        severity_mapping = {'Control': 0, 'Classical_Dengue_Fever': 1, 
                          'Dengue_Hemorrhagic_Fever': 2, 'Dengue_Shock_Syndrome': 3}
        chain_data['severity_numeric'] = chain_data['group'].map(severity_mapping)
        
        axes3[0].scatter(chain_data['severity_numeric'] + np.random.normal(0, 0.05, len(chain_data)), 
                        chain_data['expansion_gini'], alpha=0.6, label=chain, s=50)
        
        # Add trend line
        if len(chain_data) > 2:
            z = np.polyfit(chain_data['severity_numeric'], chain_data['expansion_gini'], 1)
            p = np.poly1d(z)
            x_range = np.array([0, 3])
            axes3[0].plot(x_range, p(x_range), '--', alpha=0.5)
    
    axes3[0].set_xlabel('Dengue Severity Level\n(0=Control → 3=Shock Syndrome)', fontweight='bold')
    axes3[0].set_ylabel('Gini Coefficient', fontweight='bold')
    axes3[0].set_title('Gini Coefficient vs Dengue Severity', fontweight='bold')
    axes3[0].set_xticks([0, 1, 2, 3])
    axes3[0].set_xticklabels(['Control', 'Classical', 'Hemorrhagic', 'Shock'])
    axes3[0].legend(title='Chain', bbox_to_anchor=(1.05, 1), loc='upper left')
    axes3[0].grid(True, alpha=0.3)
    
    # Violin plot
    sns.violinplot(data=gini_df, x='group', y='expansion_gini', 
                   ax=axes3[1], palette=severity_palette, order=severity_order)
    axes3[1].set_xlabel('Dengue Severity Group', fontweight='bold')
    axes3[1].set_ylabel('Gini Coefficient', fontweight='bold')
    axes3[1].set_title('Distribution by Severity Group', fontweight='bold')
    axes3[1].tick_params(axis='x', rotation=45)
    axes3[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    ResultsManager.save_figure(fig3, 'trend_analysis.png', 'plots')
    plt.close(fig3)
    
    # 4. B-cell vs T-cell comparison
    print("4. Creating B-cell vs T-cell comparison...")
    fig4, axes4 = plt.subplots(1, 2, figsize=(12, 5))
    
    # B-cells
    b_cells = gini_df[gini_df['chain'].isin(['IGH', 'IGK', 'IGL'])]
    sns.boxplot(data=b_cells, x='chain', y='expansion_gini', hue='group',
                ax=axes4[0], palette=severity_palette, order=['IGH', 'IGK', 'IGL'])
    axes4[0].set_title('B-cell Receptor Chains', fontweight='bold')
    axes4[0].set_xlabel('B-cell Chain', fontweight='bold')
    axes4[0].set_ylabel('Gini Coefficient', fontweight='bold')
    axes4[0].legend(title='Severity', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # T-cells
    t_cells = gini_df[gini_df['chain'].isin(['TRA', 'TRB'])]
    if len(t_cells) > 0:
        sns.boxplot(data=t_cells, x='chain', y='expansion_gini', hue='group',
                   ax=axes4[1], palette=severity_palette, order=['TRA', 'TRB'])
        axes4[1].set_title('T-cell Receptor Chains', fontweight='bold')
        axes4[1].set_xlabel('T-cell Chain', fontweight='bold')
        axes4[1].set_ylabel('Gini Coefficient', fontweight='bold')
        axes4[1].legend(title='Severity', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    ResultsManager.save_figure(fig4, 'bcell_tcell_comparison.png', 'plots')
    plt.close(fig4)
    
    # 5. Additional visualizations
    print("5. Creating additional visualizations...")
    
    # 5a. Distribution plot
    fig5a, ax5a = plt.subplots(figsize=(10, 6))
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]['expansion_gini'].dropna()
        if len(chain_data) > 1:
            density = gaussian_kde(chain_data)
            xs = np.linspace(0, 1, 200)
            ax5a.plot(xs, density(xs), label=chain, linewidth=2)
    
    ax5a.set_xlabel('Gini Coefficient', fontweight='bold')
    ax5a.set_ylabel('Density', fontweight='bold')
    ax5a.set_title('Distribution of Gini Coefficients by Chain', fontweight='bold')
    ax5a.legend(title='Chain')
    ax5a.grid(True, alpha=0.3)
    ax5a.set_xlim(0, 1)
    plt.tight_layout()
    ResultsManager.save_figure(fig5a, 'gini_distributions.png', 'plots')
    plt.close(fig5a)
    
    # 5b. Mean Gini by severity bar plot
    fig5b, ax5b = plt.subplots(figsize=(10, 6))
    mean_by_group = gini_df.groupby(['group', 'chain'])['expansion_gini'].mean().unstack()
    mean_by_group.plot(kind='bar', ax=ax5b, width=0.8)
    ax5b.set_xlabel('Dengue Severity Group', fontweight='bold')
    ax5b.set_ylabel('Mean Gini Coefficient', fontweight='bold')
    ax5b.set_title('Mean Gini Coefficient by Severity Group and Chain', fontweight='bold')
    ax5b.legend(title='Chain', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax5b.tick_params(axis='x', rotation=45)
    ax5b.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    ResultsManager.save_figure(fig5b, 'mean_gini_by_group.png', 'plots')
    plt.close(fig5b)
    
    print("All visualizations saved to RESULTS/plots/")

# ============================================================================
# PART 4: SUMMARY TABLES
# ============================================================================

def create_summary_tables(gini_df):
    """Create summary tables"""
    print("\n" + "="*60)
    print("PART 4: CREATING SUMMARY TABLES")
    print("="*60)
    
    # 1. Descriptive statistics by group and chain
    print("1. Creating descriptive statistics table...")
    desc_stats = gini_df.groupby(['group', 'chain'])['expansion_gini'].agg([
        'mean', 'std', 'min', 'max', 'count'
    ]).round(3)
    
    desc_stats = desc_stats.reset_index()
    ResultsManager.save_csv(desc_stats, 'descriptive_statistics.csv', 'tables')
    
    # 2. Mean Gini by severity
    print("2. Creating mean Gini by severity table...")
    mean_by_severity = gini_df.groupby('group')['expansion_gini'].agg([
        'mean', 'std', 'count'
    ]).round(3)
    
    mean_by_severity = mean_by_severity.reset_index()
    ResultsManager.save_csv(mean_by_severity, 'mean_gini_by_severity.csv', 'tables')
    
    # 3. Percentage change from control
    print("3. Creating percentage change table...")
    control_means = gini_df[gini_df['group'] == 'Control'].groupby('chain')['expansion_gini'].mean()
    
    changes = []
    for chain in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']:
        chain_data = gini_df[gini_df['chain'] == chain]
        control_mean = chain_data[chain_data['group'] == 'Control']['expansion_gini'].mean()
        
        for disease_group in ['Classical_Dengue_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']:
            disease_mean = chain_data[chain_data['group'] == disease_group]['expansion_gini'].mean()
            if not pd.isna(control_mean) and not pd.isna(disease_mean) and control_mean != 0:
                percent_change = ((disease_mean - control_mean) / control_mean) * 100
                changes.append({
                    'Chain': chain,
                    'Comparison': f'Control vs {disease_group}',
                    'Control_Mean': control_mean,
                    f'{disease_group}_Mean': disease_mean,
                    'Absolute_Change': disease_mean - control_mean,
                    'Percent_Change': percent_change
                })
    
    if changes:
        changes_df = pd.DataFrame(changes)
        ResultsManager.save_csv(changes_df, 'percentage_changes.csv', 'tables')
    
    # 4. Summary statistics by chain
    print("4. Creating chain summary table...")
    chain_summary = gini_df.groupby('chain')['expansion_gini'].agg([
        'mean', 'std', 'median', 'min', 'max', 'count'
    ]).round(3)
    
    chain_summary = chain_summary.reset_index()
    ResultsManager.save_csv(chain_summary, 'chain_summary.csv', 'tables')
    
    # 5. Statistical significance summary
    print("5. Creating statistical significance summary...")
    
    # Calculate significant differences
    significance_summary = []
    chains = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    groups = ['Control', 'Classical_Dengue_Fever', 'Dengue_Hemorrhagic_Fever', 'Dengue_Shock_Syndrome']
    
    for chain in chains:
        chain_data = gini_df[gini_df['chain'] == chain]
        
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                group1 = groups[i]
                group2 = groups[j]
                
                data1 = chain_data[chain_data['group'] == group1]['expansion_gini'].dropna()
                data2 = chain_data[chain_data['group'] == group2]['expansion_gini'].dropna()
                
                if len(data1) > 1 and len(data2) > 1:
                    _, p_value = mannwhitneyu(data1, data2)
                    
                    significance_summary.append({
                        'Chain': chain,
                        'Comparison': f'{group1} vs {group2}',
                        'p_value': p_value,
                        'Significant': 'Yes' if p_value < 0.05 else 'No',
                        'Effect_Direction': 'Increase' if data2.mean() > data1.mean() else 'Decrease' if data2.mean() < data1.mean() else 'Equal'
                    })
    
    if significance_summary:
        sig_df = pd.DataFrame(significance_summary)
        ResultsManager.save_csv(sig_df, 'significance_summary.csv', 'tables')
    
    print("All summary tables saved to RESULTS/tables/")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main analysis pipeline"""
    
    try:
        # Part 1: Calculate Gini coefficients
        gini_df = calculate_all_gini_coefficients()
        
        if gini_df is not None:
            # Part 2: Statistical analysis
            perform_statistical_analysis(gini_df)
            
            # Part 3: Visualizations
            create_main_visualizations(gini_df)
            
            # Part 4: Summary tables
            create_summary_tables(gini_df)
            
            # Final summary
            print("\n" + "="*60)
            print("ANALYSIS COMPLETED SUCCESSFULLY!")
            print("="*60)
            print(f"Analysis finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("\nAll results saved in RESULTS folder:")
            print("\nRESULTS/")
            print("├── plots/")
            print("│   ├── statistical_summary.png")
            print("│   ├── p_value_heatmap.png")
            print("│   ├── trend_analysis.png")
            print("│   ├── bcell_tcell_comparison.png")
            print("│   ├── gini_distributions.png")
            print("│   └── mean_gini_by_group.png")
            print("├── statistics/")
            print("│   ├── kruskal_wallis_results.csv")
            print("│   ├── pairwise_mann_whitney_results.csv")
            print("│   └── effect_size_analysis.csv")
            print("├── tables/")
            print("│   ├── descriptive_statistics.csv")
            print("│   ├── mean_gini_by_severity.csv")
            print("│   ├── percentage_changes.csv")
            print("│   ├── chain_summary.csv")
            print("│   └── significance_summary.csv")
            print("├── reports/")
            print("│   └── statistical_analysis_report.txt")
            print("└── raw_data/")
            print("    └── gini_coefficients.csv")
            
            # Print quick findings
            print("\nQUICK FINDINGS:")
            print("-" * 40)
            
            sig_count = 0
            for chain in ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']:
                chain_data = gini_df[gini_df['chain'] == chain]
                c_data = chain_data[chain_data['group'] == 'Control']['expansion_gini'].dropna()
                s_data = chain_data[chain_data['group'] == 'Dengue_Shock_Syndrome']['expansion_gini'].dropna()
                
                if len(c_data) > 1 and len(s_data) > 1:
                    _, p_value = mannwhitneyu(c_data, s_data)
                    if p_value < 0.05:
                        sig_count += 1
                        mean_diff = s_data.mean() - c_data.mean()
                        print(f"  {chain}: Control vs Shock: p={p_value:.4f}, Δ={mean_diff:.3f}")
            
            print(f"\nSignificant differences (Control vs Shock): {sig_count}/5 chains")
            
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please ensure all required input files are in the current directory.")
    except Exception as e:
        print(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
