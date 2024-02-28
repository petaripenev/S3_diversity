#!/usr/bin/env python3
'''Reads a vsearch cluster file with gene ids, groups them, then identifies the longest 
scaffold for each group. Outputs the longest scaffold for each group to a fasta file 
and the mapping between representative genes for clusters to the longest scaffold.'''
import csv, os
from Bio import SeqIO
from pathos.multiprocessing import ProcessingPool as Pool

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

def get_unique_names(gene_dict):
    unique_names = list()
    for group in gene_dict:
        # Iterate over each gene in the group
        for gene in gene_dict[group]:
            # Open the fasta file
            basename = "_".join(gene.split("_")[:-2])
            unique_names.append(basename)
            unique_names = list(set(unique_names))
    return unique_names

def index_fasta_file(name, location_base, scaffold_suffix):
    location = f"{location_base}{name}{scaffold_suffix}"
    # Read the fasta file
    records = SeqIO.index(location, 'fasta')
    print(f'Indexed {name}')
    return (name, records)

def get_longest_scaffold_ids(big_fasta_files_as_seqio_indexes, gene_dict):
    longest_scaffolds = {}
    for group in gene_dict:
        scaffold_lengths = {}
        for gene in gene_dict[group]:
            basename = "_".join(gene.split("_")[:-2])
            contig_name = "_".join(gene.split("_")[:-1])
            record = big_fasta_files_as_seqio_indexes[basename][contig_name]
            scaffold_lengths[record.id] = len(record.seq)
        longest_scaffold = max(scaffold_lengths, key=scaffold_lengths.get)
        longest_scaffolds[group] = longest_scaffold
    return longest_scaffolds

def main(args):

    # Load arguments
    clusters_path = args.path_to_clusters
    num_processes = args.num_processes
    output_dir = args.output_dir
    scaffold_dir = args.scaffolds_dir
    scaffolds_suffix = args.scaffolds_suffix

    # Read the cluster file and get unique names
    gene_dict = read_S3_clusters(clusters_path)
    unique_names = get_unique_names(gene_dict)

    # Add a check for available threads of the system
    available_threads = os.cpu_count()
    if num_processes > available_threads:
        print(f'Number of processes requested ({num_processes}) is greater than available threads ({available_threads}).')
        print(f'Using {available_threads} processes.')
        num_processes = available_threads

    # Index the fasta files in a multiprocessing pool
    pool = Pool(processes=num_processes)
    results = pool.map(index_fasta_file, 
                       unique_names, 
                       [scaffold_dir]*len(unique_names), 
                       [scaffolds_suffix]*len(unique_names)
            )
    pool.close()
    pool.join()
    # Store the indexes in a dictionary
    big_fasta_files_as_seqio_indexes = dict(results)

    longest_scaffolds = get_longest_scaffold_ids(big_fasta_files_as_seqio_indexes, gene_dict)
    unique_longest_scaffolds = list(set(longest_scaffolds.values()))

    # Output the mapping between cluster gene and longest scaffold
    with open(f'{output_dir}longest_scaffold_to_gene_mapping.tsv', 'w') as f:
        f.write('cluster_gene\tlongest_scaffold\n')
        for group in longest_scaffolds:
            f.write(f'{group}\t{longest_scaffolds[group]}\n')

    #Output the longest scaffold for each group to a fasta file
    with open(f'{output_dir}longest_scaffolds.fa', 'w') as f:
        for longest_id in unique_longest_scaffolds:
            basename = "_".join(longest_id.split("_")[:-1])
            record = big_fasta_files_as_seqio_indexes[basename][longest_id]
            SeqIO.write(record, f, 'fasta')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--path_to_clusters', type=str, default='results/s3/all_S3_pfam_clusters.txt', help='Path to protein cluster file')
    parser.add_argument('--num_processes', type=int, default=12, help='Number of processes to use when reading scaffold files.')
    parser.add_argument('--scaffolds_dir', type=str, default='results/scaffolds/', help='Path to directory with all scaffolds.')
    parser.add_argument('--scaffolds_suffix', type=str, default='_scaffold.fa', help='Suffix for scaffold files.')
    parser.add_argument('--output_dir', type=str, default='results/s3/', help='Output directory')

    args = parser.parse_args()
    main(args)


# for name in unique_names:
#     location = f"/groups/banfield/scratch/projects/environmental/spot/int/2023/assembly.d/pipeline/results/{name}/{name}_scaffold.fa"
#     # Read the fasta file
#     records = SeqIO.index(location, 'fasta')
#     print(f'Indexed {name}')
#     big_fasta_files_as_seqio_indexes[name] = records

