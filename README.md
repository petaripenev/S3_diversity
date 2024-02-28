# Abundance analysis of metagenomes

This snakemake workflow can be used to analyse the abundance of species groups, based on different marker proteins. 

## Inputs

The metagenomes to study are defined in the [samples.tsv](./config/samples.tsv) file and the proteins to use as markers are defined in the [proteins.tsv](./config/proteins.tsv) file.

The columns to use in the [samples.tsv](./config/samples.tsv) file are:
- **sample_name**: name of sample
- **scaffolds_path**: path to scaffolds or contigs
- **contig_to_bin**: path to the contig to bin file (generated by DAStool). Not needed for the snakemake, but can be used with the [s3_bin_completeness.py](./workflow/scripts/s3_bin_completeness.py) script.
- **forward_read_path**: path toforward reads
- **reverse_read_path**: path toreverse reads

In the same directory as **scaffolds_path** the files {scaffolds_path}.genes.fna and {scaffolds_path}.genes.faa should be present. Those should have been generated by the Snakemake assembly and annotaion pipeline.

The columns to use in the [proteins.tsv](./config/proteins.tsv) file are:

- **protein_name**: name of protein
- **hmm_path**: path to hmm file
- **filter_length_nt**: Remove protein hits with shorter length than specified.
- **cluster_threshold**: Currently not used, uses default of 99% (should use thresholds from [Matt olm et al.](https://journals.asm.org/doi/full/10.1128/mSystems.00731-19))

## Configs
In the [cluster.yaml](./config/cluster.yaml) file you can specify the way the slurm outputs are written using the `o` option. One can also set up an email to receive updates, however this also needs to be added to the command line when running snakemake.

In the [config.yaml](./config/config.yaml) file you can specify the following:
- **min_map_id**: Minimum mapping identity to use for bbmap
- **coverm_perc_id**: Minimum percentage id to use for coverage calculations

### Time and Depth
The script [s3_taxonomy_analysis.py](./workflow/scripts/s3_taxonomy_analysis.py) assumes that the sample name holds information about depth and time point. The script will use the third string in the name as time and the sixth string as depth. This can be changed with the `time_pos` and `depth_pos` parameters in the [config.yaml](./config/config.yaml) file. For now we support only a single digit.

## Test command to run

```snakemake -nrp --cluster-config config/cluster.yaml --default-resources partition=standard --cluster "sbatch -J RPclust -p {cluster.p} -o {cluster.o}" -j10  --rerun-incomplete --latency-wait 60 --cluster-cancel scancel > test_run_p.txt```

