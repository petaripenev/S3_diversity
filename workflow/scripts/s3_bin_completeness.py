#!/usr/bin/env python3
'''Analyse S3 results with the assembly and binning results'''

import sys, os
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from s3_taxonomy_analysis import read_S3_clusters, combine_s3_results
from itertools import product

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--s3_clusters', type=str, help='S3 clusters file', default='./results/S3c/all_S3c_pfam_clusters.txt')
    parser.add_argument('--s3_blast_results', type=str, help='S3 blast results file', default='./results/S3c/all_S3c_centroids_NR.tsv')
    parser.add_argument('--s3_blast_taxonomy', type=str, help='S3 blast taxonomy file', default='./results/S3c/all_S3c_centroids_NR_taxa.tsv')
    parser.add_argument('--longest_scaffold_to_gene', type=str, help='Longest scaffold to gene mapping file', default='./results/S3c/longest_scaffold_to_gene_mapping.tsv')
    parser.add_argument('--rps3_counts', type=str, help='RPS3 counts file', default='./results/maps/S3c_count_table.tsv')
    parser.add_argument('--covhist_folder', type=str, help='Folder containing coverage histograms', default='./results/maps/S3c/')
    parser.add_argument('--covstats_folder', type=str, help='Folder containing coverage stats', default='./results/maps/S3c/')
    parser.add_argument('--dastool_bins', type=str, help='DAS Tool bins folder', default='./results/dastool_bins/')
    parser.add_argument('--gene_annotations', type=str, help='Gene annotations folder', default='./results/gene_annotations/')
    parser.add_argument('--genes_of_interest', type=str, help='Genes of interest', default='xoxF,coxL,cynS,amiF,qhpA,nthA,nthB,amoA,amoB,pmoA,pmoB,nosZ,moaA,moaB,moaC,nasA,wecB,mogA,moeA')
    parser.add_argument('--color_map', type=str, help='Custom color map. List hexcolors separated by comma, guard with quotes.', default=None)
    parser.add_argument('--taxonomy_rank', type=str, help='Taxonomy rank to plot', default='phylum')
    parser.add_argument('--unassigned', help='Include unassigned taxonomies in plot', default=False, action='store_true')
    return parser.parse_args(args)

#--genes_of_interest rmlA,rmlC,wzc,wzy,pelA,epsC,betB,dnaK,otsA,ompR,proA,mscL,osmC,nasA,xoxF,coxL,moaA,cynS,nthA,pmoB --color_map "#9F2B68,#753B6B,#6B0558,#800080,#A16F89,#DE00DE,#00FFFF,#C4E2FF,#7393B3,#5D5DFC,#0000FF,#1434A4,#0E236E,#DAA520,#E49B0F,#FFBF00,#E1C16E,#E4D00A,#FDDA0D,#FAFA33"

def load_big_annotation_files(gene_annotations_folder, sample_name):
    kegg_file_name = gene_annotations_folder + sample_name + '_min1000.fa.genes.faa-vs-kegg.b6+'
    uniprot_file_name = gene_annotations_folder + sample_name + '_min1000.fa.genes.faa-vs-uniprot.b6+'
    uniref_file_name = gene_annotations_folder + sample_name + '_min1000.fa.genes.faa-vs-uni.b6+'
    annotation_files, annotations = [kegg_file_name, uniprot_file_name, uniref_file_name], list()
    for file_name in annotation_files:
        data = pd.read_csv(file_name, sep='\t', header=None)
        data.columns = ['gene', 'annotation', 'identity', 'length', 'mismatches', 'gapopenings', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'name']
        data['Contig'] = data['gene'].apply(lambda x: '_'.join(x.split('_')[:-1]))
        annotations.append(data)
    return annotations

def get_annotations_for_bin(annotations, dastool_bins_to_genes, bin_of_interest):
    for annotation in annotations:
        yield annotation[annotation['Contig'].isin(dastool_bins_to_genes[bin_of_interest])]

def get_abundances_by_gene_name(bins_of_interest, genes_of_interest, s3_results, sample_cols, tax_rank, bin_annotations):
    # Get sample abundances for each bin of interest if the gene is present
    bin_abundances = pd.DataFrame(columns = sample_cols+[tax_rank, 'gene', 'gene_presence', 'bin'])
    
    i = 0
    for bin_of_interest, gene in product(bins_of_interest, genes_of_interest):
        bin_tax = s3_results.loc[s3_results['dastool_bin'] == bin_of_interest,[tax_rank]].values[0][0]
        if bin_tax != bin_tax:
            bin_tax = 'Unassigned'
        #Search within bin_annotations name_* for the gene
        gene_presence = bin_annotations[bin_of_interest]['name_kegg'].str.contains(f'{gene}[ ;]', na=False, regex=True) | \
                        bin_annotations[bin_of_interest]['name_uniprot'].str.contains(f'{gene}[ ;]', na=False, regex=True) | \
                        bin_annotations[bin_of_interest]['name_uniref'].str.contains(f'{gene}[ ;]', na=False, regex=True)
        sample_abundance = s3_results.loc[s3_results['dastool_bin'] == bin_of_interest, sample_cols]
        if sum(gene_presence) == 0:
            continue
        bin_abundances.loc[i] = list(sample_abundance.values[0]) + [bin_tax, gene, sum(gene_presence), bin_of_interest]
        i += 1
    
    #Save the bin_abundances to file
    bin_abundances.to_csv(f'./results/bin_annotations/bin_abundances_{tax_rank}_{",".join(genes_of_interest)}.tsv', sep='\t', index=False)
    return bin_abundances

def main(args):
    cl_args = parse_args(args)
    tax_rank = cl_args.taxonomy_rank
    #Load data in pandas dataframes

    s3_clusters_dict = read_S3_clusters(cl_args.s3_clusters)
    s3_blast_results = pd.read_csv(cl_args.s3_blast_results, sep='\t', header=None)
    longest_scaffold_to_gene = pd.read_csv(cl_args.longest_scaffold_to_gene, sep='\t')
    rps3_counts = pd.read_csv(cl_args.rps3_counts, sep='\t')
    s3_blast_taxonomy = pd.read_csv(cl_args.s3_blast_taxonomy, sep='\t', header=None)
    
    #Read each file in dastool_bins folder
    dastool_genes_to_bins, dastool_bins_to_genes = dict(), dict()
    for bin_file in os.listdir(cl_args.dastool_bins):
        with open(cl_args.dastool_bins + bin_file, 'r') as f:
            #Read into dictionary with first col as key and second col as value
            for line in f:
                line = line.strip().split('\t')
                dastool_genes_to_bins[line[0]] = line[1]
                if line[1] not in dastool_bins_to_genes.keys():
                    dastool_bins_to_genes[line[1]] = list()
                dastool_bins_to_genes[line[1]].append(line[0])
    #Combine the dataframes
    s3_results, rps3_cols = combine_s3_results(s3_clusters_dict, 
                                               s3_blast_results, 
                                               longest_scaffold_to_gene, 
                                               rps3_counts, 
                                               s3_blast_taxonomy, 
                                               tax_rank)

    #Check how many of the dastool_bins keys are found in the s3_results['cluster'] list and add the dastools_bin value as new column
    dastool_bin_list = []
    for cluster in s3_results['cluster']:
        #Get any gene in the cluster that is in the dastool_bins dictionary
        bin = [dastool_genes_to_bins[gene.rsplit('_', 1)[0]] for gene in cluster if gene.rsplit('_', 1)[0] in dastool_genes_to_bins.keys()]
        if len(bin) > 1:
            #print(f'Multiple bins found for cluster: {" ".join(cluster)}')
            dastool_bin_list.append(bin[0])
        elif len(bin) == 0:
            dastool_bin_list.append('None')
        else:
            dastool_bin_list.append(bin[0])
    s3_results['dastool_bin'] = dastool_bin_list

    num_clusters_binned = sum(s3_results['dastool_bin'] != 'None')
    average_depth_binned = np.mean(s3_results[s3_results['dastool_bin'] != 'None'][list(rps3_cols.values())].sum())
    std_depth_binned = np.std(s3_results[s3_results['dastool_bin'] != 'None'][list(rps3_cols.values())].sum())
    total_depth_per_sample = np.mean(s3_results[list(rps3_cols.values())].sum())
    print(f"Number of clusters binned: {num_clusters_binned}; which is {num_clusters_binned/len(s3_results)*100:.2f}% of total clusters.")
    print(f"Percent of total depth in binned clusters: {average_depth_binned/total_depth_per_sample*100:.2f}% +/- {std_depth_binned/total_depth_per_sample*100:.2f}%")
    
    #bins_of_interest = 'WaterYear_AC3_D5_coass-16O-0-7.metabat2.6'
    #bin_of_interest = 'WaterYear_AC3_D5_coass-16O-0-7.concoct.131'
    #bin_of_interest = 'WaterYear_AC3_D5_coass-16O-0-7.maxbin2.170_sub'
    
    ### Some s3 clusters map to multiple bins, so in the Snakemake we should add the 600AA cut off after clustering ###
    #For now lets just get the unique bins
    bins_of_interest = list(s3_results[s3_results['dastool_bin'] != 'None']['dastool_bin'].unique())
    bin_annotations = dict()
    for bin_of_interest in bins_of_interest:
        sample_name = bin_of_interest.split('.')[0]
        taxonomy = s3_results.loc[s3_results['dastool_bin'] == bin_of_interest, tax_rank].values[0]
        if os.path.exists(f'./results/bin_annotations/{bin_of_interest}_{taxonomy}.tsv'):
            #print(f'Annotations for {bin_of_interest} already saved to file')
            bin_annotations[bin_of_interest] = pd.read_csv(f'./results/bin_annotations/{bin_of_interest}_{taxonomy}.tsv', sep='\t')
            continue
        if sample_name not in bin_annotations.keys():
            bin_annotations[sample_name] = load_big_annotation_files(cl_args.gene_annotations, sample_name)
        annotations = list(get_annotations_for_bin(bin_annotations[sample_name], dastool_bins_to_genes, bin_of_interest))

        #annotations = list(get_annotations_for_bin(cl_args.gene_annotations, dastool_bins_to_genes, bin_of_interest))
        #Merge the three annotations on gene, keeping only the gene, name, and annotation columns
        all_annotations = annotations[0].merge(annotations[1], on='gene', how='outer', suffixes=('_kegg', '_uniprot')).merge(annotations[2], on='gene', how='outer')
        all_annotations = all_annotations[['gene', 'name_kegg', 'annotation_kegg', 'name_uniprot', 'annotation_uniprot', 'name', 'annotation']]
        all_annotations.columns = ['name_gene', 'name_kegg', 'annotation_kegg', 'name_uniprot', 'annotation_uniprot', 'name_uniref', 'annotation_uniref']
        #save the annotations to a file
        all_annotations.to_csv(f'./results/bin_annotations/{bin_of_interest}_{taxonomy}.tsv', sep='\t')
        print(f'Annotations for {bin_of_interest} saved to file')
        bin_annotations[bin_of_interest] = all_annotations


    # Genes of interest
    # genes_of_interest = ['xoxF', 'coxL', 'cynS', 'amiF', 'qhpA', 'nthA', 'nthB', 'amoA', 'amoB', 'pmoA', 
    #                      'pmoB', 'nosZ', 'moaA', 'moaB', 'moaC', 'nasA', 'wecB', 'mogA', 'moeA']
    genes_of_interest = cl_args.genes_of_interest.split(',')
    
    sample_cols = list(sorted(rps3_cols.values(), reverse=True))

    #Create a gene to color mapping from tab20 or custom
    gene_colors = dict()
    if not cl_args.color_map:
        for i, gene in enumerate(genes_of_interest):
            gene_colors[gene] = plt.cm.tab20(i)
    else:
        for i, gene in enumerate(genes_of_interest):
            gene_colors[gene] = cl_args.color_map.split(',')[i]

    #Get the abundances of the genes of interest
    if os.path.exists(f'./results/bin_annotations/bin_abundances_{tax_rank}_{cl_args.genes_of_interest}.tsv'):
        bin_abundances = pd.read_csv(f'./results/bin_annotations/bin_abundances_{tax_rank}_{cl_args.genes_of_interest}.tsv', sep='\t')
    else:
        bin_abundances = get_abundances_by_gene_name(bins_of_interest, genes_of_interest, s3_results, sample_cols, tax_rank, bin_annotations)


    grouped_df = bin_abundances.groupby([tax_rank, 'gene']).sum().reset_index()
    if 'gene_presence' in grouped_df.columns:
        grouped_df.drop(columns=['gene_presence'], inplace=True)
    #Rename Candidatus to nothing
    grouped_df[tax_rank] = grouped_df[tax_rank].apply(lambda x: x.replace('Candidatus ', ''))
    #Remove unassigned
    if not cl_args.unassigned:
        grouped_df = grouped_df[grouped_df[tax_rank] != 'Unassigned']
    #Get number of unique phyla
    num_phyla = len(grouped_df[tax_rank].unique())


    dfs_by_sample = dict()
    for sample in sample_cols:
        dfs_by_sample[sample] = grouped_df[[sample, tax_rank, 'gene']].pivot(index=tax_rank, columns='gene', values=sample).fillna(0)

    # #Plot the abundances of the genes by phylum and sample
    # fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(15, 6))
    # hatches = {'D1': '','D3': '///', 'D5': '...'}
    # for i, sample in enumerate(sample_cols):
    #     #Create matplotlib colormap object from the gene_colors dictionary
    #     cmap = matplotlib.colors.ListedColormap([gene_colors.get(x) for x in dfs_by_sample[sample].columns])
    #     ax = dfs_by_sample[sample].plot(kind='bar', 
    #                                     stacked=True, 
    #                                     ax=axes, 
    #                                     position=i, 
    #                                     width=0.06, 
    #                                     align='edge',
    #                                     colormap=cmap,
    #                                     hatch=hatches[sample.split('_')[1]],
    #                                     legend=False
    #                                 )
    #     ax.set_xticklabels(ax.get_xticklabels(), rotation=60, ha='right')
    # l, r = plt.xlim()
    # plt.xlim(l-0.1, r+0.2)
    # markers = [plt.Line2D([0,0],[0,0],color=color, marker='s', linestyle='') for color in gene_colors.values()]
    # plt.legend(markers, gene_colors.keys(), numpoints=1, bbox_to_anchor=(1.05, 1), loc='upper left')
    # #plt.legend(sorted(genes_of_interest), bbox_to_anchor=(1.05, 1), loc='upper left')

    # plt.savefig(f'./results/{tax_rank}_gene_presence.svg', dpi=300, bbox_inches='tight')

    #Plot the abundances of the genes by phylum and sample in 3 subplots for each depth
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 12), sharex=True, sharey=True)
    ax_dict = {'D1': axes[0],'D3': axes[1], 'D5': axes[2]}
    for i, sample in enumerate(sorted(sample_cols)):
        ordered_data = dfs_by_sample[sample][genes_of_interest[::-1]]
        #Create matplotlib colormap object from the gene_colors dictionary
        cmap = matplotlib.colors.ListedColormap([gene_colors.get(x) for x in ordered_data.columns])
        ax = ordered_data.plot(kind='bar', 
                                sharex=True,
                                sharey=True,
                                stacked=True, 
                                ax=ax_dict[sample.split('_')[1]], 
                                position=i//len(ax_dict), 
                                width=0.08, 
                                align='edge',
                                colormap=cmap,
                                legend=False
                              )
    axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=60, ha='right')
    l, r = plt.xlim()
    for ax in axes:
        ax.set_xlim(l-0.1, r+0.1)
    #plt.xlim(l-0.1, r+0.1)
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='s', linestyle='') for color in gene_colors.values()]
    plt.legend(markers, gene_colors.keys(), numpoints=1, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'./results/{tax_rank}_{cl_args.genes_of_interest}_3depths.svg', dpi=300, bbox_inches='tight')
    

    #Output taxonomy, bin, and abundances of rps3_cols to csv
    # s3_results[s3_results['dastool_bin'] != 'None'].to_csv(f'results/s3/s3_{tax_rank}_bin_abundance.tsv', sep='\t', index=False, columns=[tax_rank, 'dastool_bin'] + list(rps3_cols.values()))



if __name__ == '__main__':
    main(sys.argv[1:])

