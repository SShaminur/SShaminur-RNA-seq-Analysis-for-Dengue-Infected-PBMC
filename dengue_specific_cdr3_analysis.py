# dengue_specific_cdr3_analysis.py
import pandas as pd
import numpy as np
import os
import re
from collections import Counter, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import warnings
warnings.filterwarnings('ignore')

# Your sample metadata
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

def extract_cdr3_sequences_from_file(file_path, sample_id, chain_type):
    """Extract CDR3 sequences from a single MiXCR output file"""
    
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        # Extract essential information
        sequences = []
        for idx, row in df.iterrows():
            # Extract CDR3 amino acid sequence
            cdr3_aa = str(row.get('aaSeqCDR3', ''))
            cdr3_nt = str(row.get('nSeqCDR3', ''))
            
            # Filter for productive sequences
            if (len(cdr3_aa) >= 8 and  # Minimum length
                cdr3_aa != 'nan' and 
                cdr3_aa != '' and
                not cdr3_aa.startswith('region_not_covered')):
                
                # Check if productive (starts with C, ends with W/F)
                is_productive = (cdr3_aa.startswith('C') and 
                                (cdr3_aa.endswith('W') or cdr3_aa.endswith('F')))
                
                # Extract gene information
                v_gene = str(row.get('allVHitsWithScore', ''))
                j_gene = str(row.get('allJHitsWithScore', ''))
                
                # Clean gene names
                v_gene_clean = re.sub(r'\(.*', '', v_gene) if v_gene != 'nan' else 'Unknown'
                j_gene_clean = re.sub(r'\(.*', '', j_gene) if j_gene != 'nan' else 'Unknown'
                
                sequences.append({
                    'sample_id': sample_id,
                    'chain_type': chain_type,
                    'cdr3_aa': cdr3_aa,
                    'cdr3_nt': cdr3_nt if cdr3_nt != 'nan' else '',
                    'v_gene': v_gene_clean,
                    'j_gene': j_gene_clean,
                    'read_count': float(row.get('readCount', 0)),
                    'read_fraction': float(row.get('readFraction', 0)),
                    'is_productive': is_productive,
                    'clone_id': idx
                })
        
        return sequences
        
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return []

def collect_all_cdr3_sequences(results_dir="."):
    """Collect all CDR3 sequences from all samples"""
    
    print("Collecting CDR3 sequences from all samples...")
    
    all_sequences = []
    chains_to_process = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB']
    
    for sample_id in DENGUE_METADATA['sample_id']:
        print(f"  Processing {sample_id}...")
        
        for chain in chains_to_process:
            file_path = os.path.join(results_dir, f"{sample_id}.clones_{chain}.tsv")
            
            if os.path.exists(file_path):
                sequences = extract_cdr3_sequences_from_file(file_path, sample_id, chain)
                all_sequences.extend(sequences)
    
    # Convert to DataFrame
    df = pd.DataFrame(all_sequences)
    
    # Add metadata
    metadata_df = pd.DataFrame(DENGUE_METADATA)
    df = pd.merge(df, metadata_df, on='sample_id', how='left')
    
    print(f"\nTotal sequences collected: {len(df)}")
    print(f"Unique CDR3 amino acid sequences: {df['cdr3_aa'].nunique()}")
    
    return df

def identify_dengue_specific_sequences(cdr3_df):
    """
    Identify potentially dengue-specific sequences by comparing dengue vs control groups
    """
    
    print("\nIdentifying dengue-specific CDR3 sequences...")
    
    # Separate dengue and control samples
    dengue_samples = cdr3_df[cdr3_df['group'] != 'Control']
    control_samples = cdr3_df[cdr3_df['group'] == 'Control']
    
    # Get unique CDR3 sequences in each group
    dengue_cdr3_set = set(dengue_samples['cdr3_aa'].unique())
    control_cdr3_set = set(control_samples['cdr3_aa'].unique())
    
    # Sequences found ONLY in dengue samples (not in controls)
    dengue_specific_sequences = dengue_cdr3_set - control_cdr3_set
    
    # Sequences found in BOTH dengue and controls
    common_sequences = dengue_cdr3_set.intersection(control_cdr3_set)
    
    print(f"  Unique sequences in dengue samples: {len(dengue_cdr3_set)}")
    print(f"  Unique sequences in control samples: {len(control_cdr3_set)}")
    print(f"  Dengue-specific sequences (not in controls): {len(dengue_specific_sequences)}")
    print(f"  Common sequences (in both): {len(common_sequences)}")
    
    return dengue_specific_sequences, common_sequences

def find_public_clones_across_dengue(cdr3_df, dengue_specific_sequences):
    """Find public clones shared across multiple dengue patients"""
    
    print("\nFinding public clones across dengue patients...")
    
    # Filter for dengue-specific sequences
    dengue_df = cdr3_df[cdr3_df['cdr3_aa'].isin(dengue_specific_sequences)]
    
    # Count how many patients have each CDR3 sequence
    public_clones = dengue_df.groupby('cdr3_aa').agg({
        'sample_id': 'nunique',
        'read_count': 'sum',
        'chain_type': lambda x: x.mode().iloc[0] if len(x) > 0 else 'Unknown',
        'v_gene': lambda x: x.mode().iloc[0] if len(x) > 0 else 'Unknown'
    }).reset_index()
    
    public_clones = public_clones.rename(columns={
        'sample_id': 'n_patients',
        'read_count': 'total_reads'
    })
    
    # Sort by number of patients and read count
    public_clones = public_clones.sort_values(['n_patients', 'total_reads'], ascending=False)
    
    print(f"  Public clones found: {len(public_clones)}")
    print(f"  Clones in ≥2 patients: {len(public_clones[public_clones['n_patients'] >= 2])}")
    print(f"  Clones in ≥3 patients: {len(public_clones[public_clones['n_patients'] >= 3])}")
    
    return public_clones

def analyze_sequence_enrichment(cdr3_df, dengue_specific_sequences):
    """Analyze enrichment of sequences in different dengue severity groups"""
    
    print("\nAnalyzing sequence enrichment by dengue severity...")
    
    enrichment_results = []
    
    for seq in list(dengue_specific_sequences)[:100]:  # Analyze top 100 for speed
        seq_data = cdr3_df[cdr3_df['cdr3_aa'] == seq]
        
        if len(seq_data) > 0:
            # Count occurrences in each severity group
            group_counts = seq_data['group'].value_counts()
            total_count = seq_data['read_count'].sum()
            
            # Get chain type and V gene
            chain_type = seq_data['chain_type'].mode().iloc[0] if len(seq_data['chain_type'].mode()) > 0 else 'Unknown'
            v_gene = seq_data['v_gene'].mode().iloc[0] if len(seq_data['v_gene'].mode()) > 0 else 'Unknown'
            
            enrichment_results.append({
                'cdr3_aa': seq,
                'cdr3_length': len(seq),
                'chain_type': chain_type,
                'v_gene': v_gene,
                'total_reads': total_count,
                'n_patients': seq_data['sample_id'].nunique(),
                'classical_count': group_counts.get('Classical_Dengue_Fever', 0),
                'hemorrhagic_count': group_counts.get('Dengue_Hemorrhagic_Fever', 0),
                'shock_count': group_counts.get('Dengue_Shock_Syndrome', 0)
            })
    
    enrichment_df = pd.DataFrame(enrichment_results)
    
    # Calculate which group has maximum count and enrichment score
    if len(enrichment_df) > 0:
        # Determine which group has maximum count
        def get_enriched_in(row):
            counts = {
                'Classical': row['classical_count'],
                'Hemorrhagic': row['hemorrhagic_count'],
                'Shock': row['shock_count']
            }
            return max(counts, key=counts.get)
        
        enrichment_df['enriched_in'] = enrichment_df.apply(get_enriched_in, axis=1)
        
        # Calculate severity score
        enrichment_df['severity_score'] = (
            enrichment_df['shock_count'] * 3 + 
            enrichment_df['hemorrhagic_count'] * 2 + 
            enrichment_df['classical_count'] * 1
        ) / (enrichment_df['total_reads'] + 1)
        
        enrichment_df = enrichment_df.sort_values(['severity_score', 'total_reads'], ascending=False)
    
    print(f"  Analyzed {len(enrichment_df)} dengue-specific sequences")
    
    return enrichment_df

def cluster_similar_cdr3_sequences(cdr3_df, similarity_threshold=0.7):
    """Cluster similar CDR3 sequences using sequence similarity"""
    
    print(f"\nClustering similar CDR3 sequences (threshold: {similarity_threshold})...")
    
    # Get unique sequences
    unique_seqs = cdr3_df['cdr3_aa'].unique()
    
    # Simple clustering based on string similarity
    clusters = []
    clustered_seqs = set()
    
    def sequence_similarity(seq1, seq2):
        """Calculate simple sequence similarity"""
        if len(seq1) == 0 or len(seq2) == 0:
            return 0
        
        # Use simple character matching
        matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
        max_len = max(len(seq1), len(seq2))
        return matches / max_len
    
    for i, seq1 in enumerate(unique_seqs):
        if seq1 in clustered_seqs:
            continue
        
        # Start new cluster
        cluster = [seq1]
        clustered_seqs.add(seq1)
        
        # Find similar sequences
        for j, seq2 in enumerate(unique_seqs[i+1:], i+1):
            if seq2 in clustered_seqs:
                continue
            
            if sequence_similarity(seq1, seq2) >= similarity_threshold:
                cluster.append(seq2)
                clustered_seqs.add(seq2)
        
        if len(cluster) > 1:  # Only keep clusters with multiple sequences
            clusters.append({
                'cluster_id': f"Cluster_{len(clusters)+1}",
                'representative_seq': cluster[0],
                'cluster_size': len(cluster),
                'sequences': cluster,
                'mean_length': np.mean([len(s) for s in cluster])
            })
    
    print(f"  Found {len(clusters)} clusters of similar sequences")
    print(f"  Largest cluster size: {max([c['cluster_size'] for c in clusters]) if clusters else 0}")
    
    return clusters

def search_for_dengue_epitope_motifs(cdr3_df, dengue_specific_sequences):
    """Search for common motifs in dengue-specific CDR3 sequences"""
    
    print("\nSearching for common motifs in dengue-specific CDR3s...")
    
    # Get dengue-specific sequences
    dengue_seqs = list(dengue_specific_sequences)
    
    if len(dengue_seqs) < 10:
        print("  Not enough sequences for motif analysis")
        return []
    
    # Look for common amino acid patterns
    motifs = []
    
    # Common motifs based on known antibody/TCR features
    known_motifs = {
        'GXG': 'Glycine-rich motif (common in CDR3 loops)',
        'DX[ST]': 'Aspartic acid followed by aromatic (potential antigen contact)',
        '[YW]G': 'Tryptophan/Tyrosine-Glycine (common in antigen binding)',
        'C.X{6,12}C': 'Cysteine spacing pattern',
        '[DE]X[DE]': 'Acidic patch',
        '[RK]X[RK]': 'Basic patch'
    }
    
    # Also look for position-specific patterns
    position_counts = defaultdict(lambda: defaultdict(int))
    
    for seq in dengue_seqs:
        for i, aa in enumerate(seq):
            position_counts[i][aa] += 1
    
    # Find conserved positions
    conserved_positions = []
    for pos in position_counts:
        total = sum(position_counts[pos].values())
        if total > len(dengue_seqs) * 0.3:  # Position present in >30% sequences
            most_common = max(position_counts[pos].items(), key=lambda x: x[1])
            if most_common[1] / total > 0.5:  # > 50% conservation
                conserved_positions.append({
                    'position': pos,
                    'amino_acid': most_common[0],
                    'frequency': most_common[1] / total,
                    'total_sequences': total
                })
    
    conserved_positions.sort(key=lambda x: x['frequency'], reverse=True)
    
    print(f"  Analyzed {len(dengue_seqs)} dengue-specific sequences")
    print(f"  Found {len(conserved_positions)} conserved positions")
    
    return conserved_positions

def create_cdr3_sequence_database(cdr3_df, output_dir="dengue_cdr3_database"):
    """Create a comprehensive database of dengue-specific CDR3 sequences"""
    
    print(f"\nCreating CDR3 sequence database in '{output_dir}'...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Save all CDR3 sequences
    all_sequences_file = os.path.join(output_dir, "all_cdr3_sequences.csv")
    cdr3_df.to_csv(all_sequences_file, index=False)
    print(f"  Saved all sequences to: {all_sequences_file}")
    
    # 2. Identify dengue-specific sequences
    dengue_specific, common = identify_dengue_specific_sequences(cdr3_df)
    
    # 3. Save dengue-specific sequences
    dengue_specific_df = cdr3_df[cdr3_df['cdr3_aa'].isin(dengue_specific)]
    dengue_specific_file = os.path.join(output_dir, "dengue_specific_cdr3.csv")
    dengue_specific_df.to_csv(dengue_specific_file, index=False)
    print(f"  Saved dengue-specific sequences to: {dengue_specific_file}")
    
    # 4. Find public clones
    public_clones = find_public_clones_across_dengue(cdr3_df, dengue_specific)
    public_clones_file = os.path.join(output_dir, "public_clones.csv")
    public_clones.to_csv(public_clones_file, index=False)
    print(f"  Saved public clones to: {public_clones_file}")
    
    # 5. Analyze enrichment
    enrichment_df = analyze_sequence_enrichment(cdr3_df, dengue_specific)
    enrichment_file = os.path.join(output_dir, "sequence_enrichment.csv")
    enrichment_df.to_csv(enrichment_file, index=False)
    print(f"  Saved enrichment analysis to: {enrichment_file}")
    
    # 6. Cluster similar sequences
    clusters = cluster_similar_cdr3_sequences(cdr3_df)
    clusters_file = os.path.join(output_dir, "sequence_clusters.csv")
    if clusters:
        pd.DataFrame(clusters).to_csv(clusters_file, index=False)
    else:
        pd.DataFrame().to_csv(clusters_file)
    print(f"  Saved sequence clusters to: {clusters_file}")
    
    # 7. Search for motifs
    motifs = search_for_dengue_epitope_motifs(cdr3_df, dengue_specific)
    motifs_file = os.path.join(output_dir, "conserved_motifs.csv")
    if motifs:
        pd.DataFrame(motifs).to_csv(motifs_file, index=False)
    else:
        pd.DataFrame().to_csv(motifs_file)
    print(f"  Saved conserved motifs to: {motifs_file}")
    
    # 8. Create FASTA files for downstream analysis
    create_fasta_files(cdr3_df, dengue_specific, output_dir)
    
    # 9. Generate summary report
    generate_summary_report(cdr3_df, dengue_specific, public_clones, enrichment_df, output_dir)
    
    return {
        'all_sequences': all_sequences_file,
        'dengue_specific': dengue_specific_file,
        'public_clones': public_clones_file,
        'enrichment': enrichment_file,
        'clusters': clusters_file,
        'motifs': motifs_file
    }

def create_fasta_files(cdr3_df, dengue_specific, output_dir):
    """Create FASTA files for BLAST analysis"""
    
    print("\nCreating FASTA files for sequence analysis...")
    
    # FASTA for all productive sequences
    all_fasta = []
    for idx, row in cdr3_df.iterrows():
        if row['is_productive']:
            record = SeqRecord(
                Seq(row['cdr3_aa']),
                id=f"{row['sample_id']}_{row['chain_type']}_{idx}",
                description=f"sample={row['sample_id']} chain={row['chain_type']} group={row['group']} reads={row['read_count']}"
            )
            all_fasta.append(record)
    
    all_fasta_file = os.path.join(output_dir, "all_productive_cdr3.fasta")
    if all_fasta:
        SeqIO.write(all_fasta, all_fasta_file, "fasta")
        print(f"  Created FASTA file: {all_fasta_file}")
    else:
        print(f"  No productive sequences found for FASTA")
    
    # FASTA for dengue-specific sequences
    dengue_fasta = []
    dengue_df = cdr3_df[cdr3_df['cdr3_aa'].isin(dengue_specific)]
    
    for idx, row in dengue_df.iterrows():
        if row['is_productive']:
            record = SeqRecord(
                Seq(row['cdr3_aa']),
                id=f"DENGUE_{row['sample_id']}_{row['chain_type']}_{idx}",
                description=f"dengue_specific sample={row['sample_id']} chain={row['chain_type']} group={row['group']} reads={row['read_count']}"
            )
            dengue_fasta.append(record)
    
    dengue_fasta_file = os.path.join(output_dir, "dengue_specific_cdr3.fasta")
    if dengue_fasta:
        SeqIO.write(dengue_fasta, dengue_fasta_file, "fasta")
        print(f"  Created FASTA file: {dengue_fasta_file}")
    else:
        print(f"  No dengue-specific sequences found for FASTA")
    
    # FASTA for nucleotide sequences
    nt_fasta = []
    for idx, row in cdr3_df.iterrows():
        if row['is_productive'] and row['cdr3_nt'] and len(row['cdr3_nt']) > 0:
            record = SeqRecord(
                Seq(row['cdr3_nt']),
                id=f"NT_{row['sample_id']}_{row['chain_type']}_{idx}",
                description=f"nucleotide sample={row['sample_id']} chain={row['chain_type']} group={row['group']}"
            )
            nt_fasta.append(record)
    
    nt_fasta_file = os.path.join(output_dir, "cdr3_nucleotide.fasta")
    if nt_fasta:
        SeqIO.write(nt_fasta, nt_fasta_file, "fasta")
        print(f"  Created FASTA file: {nt_fasta_file}")
    else:
        print(f"  No nucleotide sequences found for FASTA")

def generate_summary_report(cdr3_df, dengue_specific, public_clones, enrichment_df, output_dir):
    """Generate a comprehensive summary report"""
    
    print("\nGenerating summary report...")
    
    report_file = os.path.join(output_dir, "dengue_cdr3_analysis_report.txt")
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("DENGUE-SPECIFIC BCR/TCR CDR3 SEQUENCE ANALYSIS REPORT\n")
        f.write("="*80 + "\n\n")
        
        # Overview
        f.write("1. OVERVIEW\n")
        f.write("-"*40 + "\n")
        f.write(f"Total samples analyzed: {len(DENGUE_METADATA['sample_id'])}\n")
        f.write(f"  - Control: 6 samples\n")
        f.write(f"  - Classical Dengue Fever: 6 samples\n")
        f.write(f"  - Dengue Hemorrhagic Fever: 6 samples\n")
        f.write(f"  - Dengue Shock Syndrome: 6 samples\n")
        f.write(f"\nTotal CDR3 sequences collected: {len(cdr3_df)}\n")
        f.write(f"Unique CDR3 amino acid sequences: {cdr3_df['cdr3_aa'].nunique()}\n")
        f.write(f"Productive sequences: {cdr3_df['is_productive'].sum()}\n")
        
        # Chain distribution
        f.write("\n2. CHAIN DISTRIBUTION\n")
        f.write("-"*40 + "\n")
        chain_counts = cdr3_df['chain_type'].value_counts()
        for chain, count in chain_counts.items():
            f.write(f"  {chain}: {count} sequences\n")
        
        # Dengue-specific sequences
        f.write("\n3. DENGUE-SPECIFIC SEQUENCES\n")
        f.write("-"*40 + "\n")
        f.write(f"Sequences found ONLY in dengue patients (not in controls): {len(dengue_specific)}\n")
        
        # Public clones
        f.write("\n4. PUBLIC CLONES (Shared across patients)\n")
        f.write("-"*40 + "\n")
        f.write(f"Total public clones: {len(public_clones)}\n")
        if len(public_clones) > 0:
            f.write(f"Clones shared by ≥2 patients: {len(public_clones[public_clones['n_patients'] >= 2])}\n")
            f.write(f"Clones shared by ≥3 patients: {len(public_clones[public_clones['n_patients'] >= 3])}\n")
            
            # Top public clones
            f.write("\nTop 10 public clones:\n")
            for i, (_, row) in enumerate(public_clones.head(10).iterrows(), 1):
                seq_display = row['cdr3_aa'][:20] + "..." if len(row['cdr3_aa']) > 20 else row['cdr3_aa']
                f.write(f"  {i}. {seq_display} (patients: {row['n_patients']}, reads: {row['total_reads']}, chain: {row['chain_type']})\n")
        else:
            f.write("No public clones found.\n")
        
        # Severity enrichment
        f.write("\n5. SEVERITY ENRICHMENT\n")
        f.write("-"*40 + "\n")
        if len(enrichment_df) > 0:
            # Get top sequences for each severity
            if 'enriched_in' in enrichment_df.columns:
                top_severe = enrichment_df[enrichment_df['enriched_in'] == 'Shock'].head(5) if len(enrichment_df[enrichment_df['enriched_in'] == 'Shock']) > 0 else pd.DataFrame()
                top_moderate = enrichment_df[enrichment_df['enriched_in'] == 'Hemorrhagic'].head(5) if len(enrichment_df[enrichment_df['enriched_in'] == 'Hemorrhagic']) > 0 else pd.DataFrame()
                top_mild = enrichment_df[enrichment_df['enriched_in'] == 'Classical'].head(5) if len(enrichment_df[enrichment_df['enriched_in'] == 'Classical']) > 0 else pd.DataFrame()
                
                if len(top_severe) > 0:
                    f.write("Top sequences enriched in Shock Syndrome:\n")
                    for i, (_, row) in enumerate(top_severe.iterrows(), 1):
                        seq_display = row['cdr3_aa'][:15] + "..." if len(row['cdr3_aa']) > 15 else row['cdr3_aa']
                        f.write(f"  {i}. {seq_display} (reads: {row['total_reads']}, patients: {row['n_patients']})\n")
        
        # Recommendations
        f.write("\n6. RECOMMENDATIONS FOR FURTHER ANALYSIS\n")
        f.write("-"*40 + "\n")
        f.write("1. BLAST search against viral antigen databases\n")
        f.write("2. Structural modeling of top public clones\n")
        f.write("3. Tetramer staining for validation\n")
        f.write("4. Functional assays for neutralizing antibodies\n")
        f.write("5. Longitudinal tracking in convalescent patients\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("ANALYSIS COMPLETED\n")
        f.write("="*80 + "\n")
    
    print(f"  Summary report saved to: {report_file}")

def visualize_cdr3_analysis(cdr3_df, public_clones, output_dir):
    """Create visualizations of CDR3 analysis"""
    
    print("\nCreating visualizations...")
    
    try:
        # Set style
        plt.style.use('seaborn-whitegrid')
        
        # Create figure
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        # Colors for groups
        colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']
        
        # Plot 1: CDR3 length distribution for Control
        ax1 = axes[0]
        control_data = cdr3_df[cdr3_df['group'] == 'Control']
        if len(control_data) > 0:
            lengths = control_data['cdr3_aa'].apply(len)
            ax1.hist(lengths, bins=20, alpha=0.7, color=colors[0], edgecolor='black')
            ax1.set_title('Control\nCDR3 Length Distribution', fontweight='bold')
            ax1.set_xlabel('CDR3 Length (amino acids)')
            ax1.set_ylabel('Frequency')
            ax1.grid(True, alpha=0.3)
        
        # Plot 2: CDR3 length distribution for Classical Dengue
        ax2 = axes[1]
        classical_data = cdr3_df[cdr3_df['group'] == 'Classical_Dengue_Fever']
        if len(classical_data) > 0:
            lengths = classical_data['cdr3_aa'].apply(len)
            ax2.hist(lengths, bins=20, alpha=0.7, color=colors[1], edgecolor='black')
            ax2.set_title('Classical Dengue Fever\nCDR3 Length Distribution', fontweight='bold')
            ax2.set_xlabel('CDR3 Length (amino acids)')
            ax2.set_ylabel('Frequency')
            ax2.grid(True, alpha=0.3)
        
        # Plot 3: CDR3 length distribution for Hemorrhagic Fever
        ax3 = axes[2]
        hemorrhagic_data = cdr3_df[cdr3_df['group'] == 'Dengue_Hemorrhagic_Fever']
        if len(hemorrhagic_data) > 0:
            lengths = hemorrhagic_data['cdr3_aa'].apply(len)
            ax3.hist(lengths, bins=20, alpha=0.7, color=colors[2], edgecolor='black')
            ax3.set_title('Dengue Hemorrhagic Fever\nCDR3 Length Distribution', fontweight='bold')
            ax3.set_xlabel('CDR3 Length (amino acids)')
            ax3.set_ylabel('Frequency')
            ax3.grid(True, alpha=0.3)
        
        # Plot 4: CDR3 length distribution for Shock Syndrome
        ax4 = axes[3]
        shock_data = cdr3_df[cdr3_df['group'] == 'Dengue_Shock_Syndrome']
        if len(shock_data) > 0:
            lengths = shock_data['cdr3_aa'].apply(len)
            ax4.hist(lengths, bins=20, alpha=0.7, color=colors[3], edgecolor='black')
            ax4.set_title('Dengue Shock Syndrome\nCDR3 Length Distribution', fontweight='bold')
            ax4.set_xlabel('CDR3 Length (amino acids)')
            ax4.set_ylabel('Frequency')
            ax4.grid(True, alpha=0.3)
        
        # Plot 5: Public clone distribution
        ax5 = axes[4]
        if len(public_clones) > 0:
            patient_counts = public_clones['n_patients'].value_counts().sort_index()
            ax5.bar(patient_counts.index, patient_counts.values, color='steelblue', alpha=0.7)
            ax5.set_title('Public Clone Distribution', fontweight='bold')
            ax5.set_xlabel('Number of Patients Sharing Clone')
            ax5.set_ylabel('Number of Clones')
            ax5.grid(True, alpha=0.3, axis='y')
        
        # Plot 6: Chain usage in dengue-specific sequences
        ax6 = axes[5]
        dengue_specific = cdr3_df[~cdr3_df['cdr3_aa'].isin(
            cdr3_df[cdr3_df['group'] == 'Control']['cdr3_aa'].unique()
        )]
        if len(dengue_specific) > 0:
            chain_usage = dengue_specific['chain_type'].value_counts()
            if len(chain_usage) > 0:
                ax6.pie(chain_usage.values, labels=chain_usage.index, autopct='%1.1f%%')
                ax6.set_title('Chain Usage in Dengue-Specific Sequences', fontweight='bold')
        
        # Adjust layout and save
        plt.tight_layout()
        output_path = os.path.join(output_dir, 'cdr3_analysis_summary.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Visualizations saved to: {output_path}")
        
    except Exception as e:
        print(f"  Error creating visualizations: {e}")

# Main execution
def main():
    """Main analysis pipeline"""
    
    print("="*80)
    print("DENGUE-SPECIFIC BCR/TCR CDR3 SEQUENCE ANALYSIS")
    print("="*80)
    
    # Set output directory
    output_dir = "dengue_cdr3_database"
    
    try:
        # Step 1: Collect all CDR3 sequences
        print("\nSTEP 1: Collecting CDR3 sequences from all samples...")
        cdr3_df = collect_all_cdr3_sequences()
        
        if len(cdr3_df) == 0:
            print("\nERROR: No CDR3 sequences found!")
            print("Please check that your MiXCR output files are in the correct directory.")
            print("Expected file format: sample_id.clones_CHAIN.tsv (e.g., gA1_S13.clones_IGH.tsv)")
            return
        
        # Step 2: Create comprehensive database
        print("\nSTEP 2: Creating CDR3 sequence database...")
        database_files = create_cdr3_sequence_database(cdr3_df, output_dir)
        
        # Step 3: Load public clones for visualization
        print("\nSTEP 3: Loading results for visualization...")
        try:
            public_clones = pd.read_csv(database_files['public_clones'])
        except:
            public_clones = pd.DataFrame()
        
        # Step 4: Create visualizations
        print("\nSTEP 4: Creating visualizations...")
        visualize_cdr3_analysis(cdr3_df, public_clones, output_dir)
        
        # Step 5: Final summary
        print("\n" + "="*80)
        print("ANALYSIS COMPLETED SUCCESSFULLY!")
        print("="*80)
        
        print(f"\nOutput files saved in directory: '{output_dir}/'")
        print("\nKey files generated:")
        print("1. all_cdr3_sequences.csv - All CDR3 sequences with metadata")
        print("2. dengue_specific_cdr3.csv - Sequences found only in dengue patients")
        print("3. public_clones.csv - Clones shared across multiple patients")
        print("4. sequence_enrichment.csv - Enrichment analysis by severity")
        print("5. sequence_clusters.csv - Clusters of similar sequences")
        print("6. conserved_motifs.csv - Conserved amino acid patterns")
        print("7. *.fasta - FASTA files for BLAST analysis")
        print("8. dengue_cdr3_analysis_report.txt - Comprehensive summary report")
        print("9. cdr3_analysis_summary.png - Visual summary of findings")
        
        print("\nNext steps:")
        print("1. Run BLAST on FASTA files against viral databases")
        print("2. Validate top public clones with experimental assays")
        print("3. Correlate with clinical severity markers")
        print("4. Compare with known dengue antibody sequences")
        
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()