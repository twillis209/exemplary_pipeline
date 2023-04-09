from itertools import chain

trait_pairs = list(chain(*[[f"{traits[i]}_and_{traits[j]}" for j in range(i+1,len(traits))] for i in range(len(traits))]))

rule run_gwas_pipeline_on_traits:
    input:
        [f"results/processed_gwas/{x}.tsv.gz" for x in traits]

rule run_sumher_on_traits:
    input:
        [f"results/ldak/ldak-thin/{x}/{{variant_set}}/{{variant_type}}/sumher.cors.full" for x in trait_pairs]
    output:
        "results/ldak/ldak-thin/combined/{variant_set}/{variant_type}/compiled_results.tsv"
    run:
        compile_sumher_files(input, output[0])
