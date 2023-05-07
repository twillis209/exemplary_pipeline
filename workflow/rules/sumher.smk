import re

rule subset_variants:
    input:
        "resources/1000g/hg38/eur/qc/{variant_set}/{variant_type}/{chr}.bed",
        "resources/1000g/hg38/eur/qc/{variant_set}/{variant_type}/{chr}.bim",
        "resources/1000g/hg38/eur/qc/{variant_set}/{variant_type}/{chr}.fam",
        range_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/matching_ids/{chr}.txt"
    output:
        temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.bed"),
        temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.bim"),
        temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.fam")
    params:
        input_stem = "resources/1000g/hg38/eur/qc/{variant_set}/{variant_type}/{chr}",
        output_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}"
    threads: 4
    resources:
        mem_mb=get_mem_mb
    group: "sumher"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.input_stem} --snps-only --make-bed --extract {input.range_file} --silent --out {params.output_stem}"

rule thin_predictors:
    input:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.bed",
        "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.bim",
        "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.fam"
    output:
        thin_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/thin.in"),
        weights_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/weights.thin")
    log:
        log_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/thin.log"
    params:
        input_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}",
        output_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/thin"
    group: "sumher"
    shell:
        """
        ldak --thin {params.output_stem} --bfile {params.input_stem} --window-prune .98 --window-kb 100 > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

rule calculate_ldak_thin_taggings_for_chromosome:
    input:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.bed",
        "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.bim",
        "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}.fam",
        weights_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/weights.thin"
    output:
        tagging_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/{chr}.tagging")
    log:
        log_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/{chr}.tagging.log"
    params:
        input_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}",
        output_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/{chr}/{chr}"
    group: "sumher"
    shell:
        "ldak --calc-tagging {params.output_stem} --bfile {params.input_stem} --weights {input.weights_file} --chr {wildcards.chr} --window-kb 1000 --power -.25 > {log.log_file}"

rule join_ldak_thin_taggings:
    input:
        [f"results/merged_gwas/{{trait_A}}_and_{{trait_B}}/{{variant_set}}/{{variant_type}}/ldak/chr{x}/chr{x}.tagging" for x in range(1, 23)]
    output:
        wg_tagging_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/whole_genome.tagging"),
        chrom_taggings_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/taggings.txt")
    log:
        log_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/whole_genome.tagging.log"
    params:
        output_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/whole_genome"
    group: "sumher"
    shell:
        """
        for x in {input}; do
            echo $x >> {output.chrom_taggings_file}
        done;

        ldak --join-tagging {params.output_stem} --taglist {output.chrom_taggings_file} > {log.log_file}
        """

rule process_sum_stats:
    input:
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/merged.tsv.gz",
        metadata_file = "resources/gwas/metadata/metadata.tsv"
    output:
        gwas_file_A = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{trait_A}.assoc"),
        gwas_file_B = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{trait_B}.assoc")
    params:
        trait_A = lambda wildcards: wildcards.trait_A,
        trait_B = lambda wildcards: wildcards.trait_B,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_a_col = 'BETA.A',
        beta_b_col = 'BETA.B',
        se_a_col = 'SE.A',
        se_b_col = 'SE.B'
    threads: 8
    resources:
        runtime = 20,
        mem_mb = get_mem_mb,
        tmpdir = "tmp"
    group: "sumher"
    conda: "../envs/exemplary_pipeline.yaml"
    script:
        "../scripts/process_sum_stats.R"

rule estimate_rg_with_ldak_thin:
    input:
        wg_tagging_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/ldak/whole_genome.tagging",
        gwas_file_A = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{trait_A}.assoc",
        gwas_file_B = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{trait_B}.assoc"
    output:
        cors_full_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/sumher.cors.full"
#        progress_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/sumher.progress",
#        cors_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/sumher.cors",
#        labels_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/sumher.labels",
#        overlap_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/sumher.overlap"
    log:
        log_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/sumher.log"
    params:
        output_stem = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/sumher"
    resources:
        runtime = 10
    group: "sumher"
    shell:
        """
        ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.gwas_file_A} --summary2 {input.gwas_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """
