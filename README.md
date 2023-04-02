# Ancillary `snakemake` pipeline for the manuscript 'Accurate detection of shared genetic architecture...'

'...from GWAS summary statistics in the small-sample context', [link here](https://doi.org/10.1101/2022.10.13.512103). The main pipeline for the simulated and 2018 UKB GWAS data set analyses can be found [here](github.com/gps_paper_pipeline). 

NB: Despite the name, this is far from being an 'exemplary' `snakemake` pipeline. Instead the name comes from the status of the data sets it was created to process which, by virtue of their sample size, could produce 'exemplary' genetic correlation estimates, at least when compared to the 2018 UKB GWAS which had relatively few cases for quite a few phenotypes.

## Running the pipeline

The target `run_sumher_on_traits` will compute and compile the estimates. Depending on the `snp_set` specified, `all` or `sans_mhc`, these estimates will be computed using data sets with or without the MHC. Running from the root directory of the project, you can invoke it as follows to produce the with-MHC estimates:

```
snakemake --profile . results/ldak/ldak-thin/combined/all/compiled_results.tsv
```

This assumes `config.yaml` specifies the right profile for you, of course. You can pass `--cluster ''` to override the cluster settings and run locally.



## `conda` environment

See `environment.yaml` for the `conda` environment (name `exemplary_pipeline`) I've had active whilst running this pipeline. I've used `snakemake` version `7.25.0`, in particular some of the more up-to-date features relating to the allocation of resources to groups, so I do recommend using `>= v7.25.0`. 

