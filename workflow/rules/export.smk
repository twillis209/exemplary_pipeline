traits = ["hyperpara", "ra", "hypothy", "asthma-ex", "sle", "crohns", "derm-ecz", "t1d-cooper", "dyschr-vit", "uc-delange", "hyperchol-ex"]

rule run_sumher_on_traits:
    input:
        [f"results/ldak/ldak-thin/{{trait}}_and_{x}/{{snp_set}}/sumher.cors.full" for x in traits]
    output:
        "results/ldak/ldak-thin/combined/{snp_set}/compiled_{trait}_results.tsv"
    run:
        compile_sumher_files(input, output[0])

