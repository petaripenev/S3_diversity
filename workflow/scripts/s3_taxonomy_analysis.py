#!/usr/bin/env python3
'''Analyse S3 results with the assembly and binning results'''

import sys, csv
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--s3_clusters', type=str, help='S3 clusters file', default='./results/S3c/all_S3c_pfam_clusters.txt')
    parser.add_argument('--s3_blast_results', type=str, help='S3 blast results file', default='./results/S3c/all_S3c_centroids_NR.tsv')
    parser.add_argument('--s3_blast_taxonomy', type=str, help='S3 blast taxonomy file', default='./results/S3c/all_S3c_centroids_NR_taxa.tsv')
    parser.add_argument('--longest_scaffold_to_gene', type=str, help='Longest scaffold to gene mapping file', default='./results/S3c/longest_scaffold_to_gene_mapping.tsv')
    parser.add_argument('--rps3_counts', type=str, help='RPS3 counts file', default='./results/maps/S3c_count_table.tsv')
    parser.add_argument('--taxonomy_rank', type=str, help='Taxonomy rank to plot', default='phylum')
    parser.add_argument('--unassigned', help='Include unassigned taxonomies in plot', default=False, action='store_true')
    parser.add_argument('--output', type=str, help='Output file name', default='s3_taxonomy_barplot.png')
    parser.add_argument('--axis_protein', type=str, help='Name of protein to plot on x-axis', default='S3c')
    return parser.parse_args(args)

def read_S3_clusters(path_to_clusters):
    gene_dict = dict()
    with open(path_to_clusters, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == 'C':
                continue
            if row[0] == 'S':
                gene_dict[row[8]] = [row[8]]
            if row[0] == 'H':
                gene_dict[row[9]].append(row[8])
    return gene_dict

def plot_stacked_bars(unassigned_rank, tax_rank, ggkbase_phylum_colors, grouped_df, prot_name, color_cols):
    if unassigned_rank:
        # Rename NaN to 'Unclassified'
        cols_grouped_df = pd.Series(grouped_df.columns)
        cols_grouped_df = cols_grouped_df.fillna('Unclassified')
        color_cols = [x if x == x else 'Unclassified' for x in color_cols]
        grouped_df.columns = cols_grouped_df
        ggkbase_phylum_colors = pd.concat([ggkbase_phylum_colors, pd.DataFrame({'Phylum': 'Unclassified', 'Color': '#000000'}, index=[0])])

    custom_cmap = list(cm.get_cmap('tab20').colors)
    if tax_rank == 'phylum':
        custom_cmap = list()
        for ph in color_cols:
            custom_cmap.append(ggkbase_phylum_colors[ggkbase_phylum_colors['Phylum'] == ph]['Color'].tolist()[0])

    #Set a nicer style
    sns.set_style("whitegrid")

    df_plot = grouped_df.plot(kind='bar', stacked=True, figsize=(10, 6), color=custom_cmap, edgecolor='black')
    plt.xlabel(tax_rank)
    plt.ylabel('Relative Abundance')
    plt.title(f'Relative {prot_name} Abundance of {tax_rank} per Sample')

    #Place the legend outside to the right of the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    #Rotate x-axis labels
    plt.xticks(rotation=45)
    plt.tight_layout()
    return df_plot

def combine_s3_results(s3_clusters_dict, s3_blast_results, longest_scaffold_to_gene, rps3_counts, s3_blast_taxonomy, tax_rank):

    s3_blast_results.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                'sstart', 'send', 'evalue', 'bitscore', 'taxid']

    
    s3_blast_taxonomy.columns = ['qseqid', 'taxid', 'taxid2', 'taxon_name', 'taxon_ranks']
    s3_blast_taxonomy = s3_blast_taxonomy.drop(columns=['taxid', 'taxid2'])

    #Split the rank by ; get the index of the rank called 'phylum'
    s3_blast_taxonomy['taxon_ranks'] = s3_blast_taxonomy['taxon_ranks'].str.split(';')
    s3_blast_taxonomy['taxons'] = s3_blast_taxonomy['taxon_name'].str.split(';')
    position_of_tax_rank = s3_blast_taxonomy['taxon_ranks'].apply(lambda x: x.index(tax_rank) if tax_rank in x else np.nan)
    #Get the phylum from the taxons column using the index of the phylum in the ranks column
    s3_blast_taxonomy[tax_rank] = [x[int(i)] if not np.isnan(i) else np.nan for x, i in zip(s3_blast_taxonomy['taxons'], position_of_tax_rank)]

    #Rename columns in rps3_counts
    rps3_cols = dict()
    for col in rps3_counts.columns:
        if col == 'Contig':
            continue
        rps3_cols[col] = '_'.join(col.split('_')[2:4])
        rps3_counts.rename(columns=rps3_cols, inplace=True)

    #Remove duplicates from s3_blast_results on first column, pick the one with the highest bitscore (12th column)
    s3_blast_results = s3_blast_results.sort_values(by=['qseqid', 'bitscore'], ascending=False).drop_duplicates(subset=['qseqid'], keep='first')
    #Remove duplicates from s3_blast_taxonomy on first column
    s3_blast_taxonomy = s3_blast_taxonomy.drop_duplicates(subset=['qseqid'], keep='first')

    #Map taxon_name to s3_blast_results from s3_blast_taxonomy using qseqid as key
    s3_blast_results[tax_rank] = s3_blast_results['qseqid'].map(s3_blast_taxonomy.set_index('qseqid')[tax_rank])

    #Check that the sorted rps3_counts Contig column is the same as the sorted longest_scaffold_to_gene longest_scaffold column
    if set(longest_scaffold_to_gene.sort_values(by=['longest_scaffold'])['longest_scaffold']) != set(rps3_counts.sort_values(by=['Contig'])['Contig']):
        print(f'ERROR: Columns Contig in rps3_counts and longest_scaffold in longest_scaffold_to_gene should match!')
        sys.exit(1)

######### CHECK THE ISSUES HERE! #########
    #Merge rps3_counts with longest_scaffold_to_gene
    rps3_counts['cluster_gene'] = rps3_counts['Contig'].map(longest_scaffold_to_gene.set_index('longest_scaffold')['cluster_gene'].to_dict())
    # if set(s3_blast_results['qseqid']) != set(rps3_counts['cluster_gene']):
    #     print(f'ERROR: Columns Contig in rps3_counts and longest_scaffold in longest_scaffold_to_gene should match!')
    #     sys.exit(1)

    #Merge rps3_counts with s3_blast_results with pd.merge on qseqid and cluster_gene
    s3_results = pd.merge(s3_blast_results, rps3_counts, how='inner', left_on='qseqid', right_on='cluster_gene')
    if sum(s3_results['cluster_gene'] != s3_results['qseqid']) != 0:
        print(f'ERROR: Columns cluster_gene and qseqid should match!')
        sys.exit(1)

    #Add a column with the clusters from s3_clusters_dict to the df
    s3_results['cluster'] = s3_results['qseqid'].map(s3_clusters_dict)

    return s3_results, rps3_cols

def main(args):
    cl_args = parse_args(args)
    tax_rank = cl_args.taxonomy_rank
    #Load data in pandas dataframes

    s3_clusters_dict = read_S3_clusters(cl_args.s3_clusters)
    s3_blast_results = pd.read_csv(cl_args.s3_blast_results, sep='\t', header=None)
    longest_scaffold_to_gene = pd.read_csv(cl_args.longest_scaffold_to_gene, sep='\t')
    rps3_counts = pd.read_csv(cl_args.rps3_counts, sep='\t')
    s3_blast_taxonomy = pd.read_csv(cl_args.s3_blast_taxonomy, sep='\t', header=None)
    ggkbase_phylum_colors = pd.read_csv('./config/ggkbase_color_scheme_phylum.csv', sep=',')

    #Combine the dataframes
    s3_results, rps3_cols = combine_s3_results(s3_clusters_dict, 
                                               s3_blast_results, 
                                               longest_scaffold_to_gene, 
                                               rps3_counts, 
                                               s3_blast_taxonomy, 
                                               tax_rank)

    # Group the s3_results dataframe by the taxonomic rank and the biological sample columns
    grouped_df = s3_results.groupby(tax_rank).sum()
    unassigned_taxa = 'noUnclassified'
    if cl_args.unassigned:
        unassigned_taxa = 'yesUnclassified'
        grouped_df = s3_results.groupby(tax_rank, dropna=False).sum()
    grouped_df = grouped_df[rps3_cols.values()]
    # Calculate the total abundance for each biological sample
    total_abundance = grouped_df.sum(axis=0)

    # Divide the abundance of each taxonomic rank by the total abundance for each biological sample
    grouped_df = grouped_df.div(total_abundance, axis=1)

    # Get the top 20 phyla by abundance
    top_phyla = grouped_df.sum(axis=1).sort_values(ascending=False).head(20).index.tolist()

    # Subset the dataframe to only include the top 20 phyla
    grouped_df = grouped_df.loc[top_phyla]
    
    #Transpose the dataframe
    grouped_df = grouped_df.T
    
    #Split sample names by _ and get the first part of the split assign as new column
    color_cols = grouped_df.columns.tolist()
    grouped_df['Time'] = grouped_df.index.str.split('').str[3]
    grouped_df['Depth'] = grouped_df.index.str.split('').str[6]

    #Sort by depth then time
    grouped_df = grouped_df.sort_values(by=['Depth', 'Time'])

    # Plot a stacked column plot of the taxonomy for each biological sample values of rps3_cols
    df_plot = plot_stacked_bars(cl_args.unassigned, tax_rank, ggkbase_phylum_colors, grouped_df, cl_args.axis_protein, color_cols)
    plt.savefig(cl_args.output, dpi=300)


if __name__ == '__main__':
    main(sys.argv[1:])

    # Split df into two dfs, one for each depth
    # grouped_df_0 = grouped_df[grouped_df['Depth'] == '3']
    # grouped_df_1 = grouped_df[grouped_df['Depth'] == '5']
    # Plot the stacked grouped barplot
    # plot_clustered_stacked([grouped_df_0, grouped_df_1], custom_cmap,["20-30cm", "50-80cm"])

    # def plot_clustered_stacked(dfall, cmap, labels=None, title="multiple stacked bar plot",  H="/", **kwargs):
    # """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
    # labels is a list of the names of the dataframe, used for the legend
    # title is a string for the title of the plot
    # H is the hatch used for identification of the different dataframe"""

    # n_df = len(dfall)
    # n_col = len(dfall[0].columns) 
    # n_ind = len(dfall[0].index)
    # axe = plt.subplot(111)

    # for df in dfall : # for each data frame
    #     axe = df.plot(kind="bar",
    #                   color=cmap,
    #                   edgecolor='black',
    #                   stacked=True,
    #                   ax=axe,
    #                   legend=False,
    #                   **kwargs)  # make bar plots

    # h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    # for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
    #     for j, pa in enumerate(h[i:i+n_col]):
    #         for rect in pa.patches: # for each index
    #             rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
    #             rect.set_hatch(H * int(i / n_col)) #edited part     
    #             rect.set_width(1 / float(n_df + 1))

    # axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    # axe.set_xticklabels(df.index, rotation = 0)
    # axe.set_title(title)

    # # Add invisible data to add another legend
    # n=[]        
    # for i in range(n_df):
    #     n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    # l1 = axe.legend(h[:n_col], l[:n_col], loc=[1.01, 0])
    # if labels is not None:
    #     l2 = plt.legend(n, labels, loc=[1.01, -0.2]) 
    # axe.add_artist(l1)
    # return axe
