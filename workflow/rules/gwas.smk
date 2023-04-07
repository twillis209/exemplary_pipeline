def get_url(w):
    url = metadata_daf.loc[metadata_daf['abbrv'] == w.trait, 'URL'].values[0]

    return url

rule download_gwas:
    output:
        "resources/gwas/{trait}.tsv.gz"
    params:
        url = get_url
    resources:
        runtime = 8
    group: "gwas"
    shell:
        """
        if [ {params.url} ]; then
            wget -O {output} {params.url}
        else
            exit -1 
        fi;
        """

rule miscellaneous_preprocessing:
    input:
        "resources/gwas/{trait}.tsv.gz"
    output:
        temp("results/processed_gwas/{trait}_pre_pipeline.tsv.gz")
    params:
        # NB: mostly from the Pan-UKB data sets
        columns_to_drop = ["af_cases_meta_hq", "af_controls_meta_hq", "beta_meta_hq", "se_meta_hq", "pval_meta_hq", "pval_heterogeneity_hq", "af_cases_meta", "af_controls_meta", "beta_meta", "se_meta", "pval_meta", "pval_heterogeneity", "af_cases_AFR", "af_cases_AMR", "af_cases_CSA", "af_cases_EAS", "af_cases_EUR", "af_cases_MID", "af_controls_AFR", "af_controls_AMR", "af_controls_CSA", "af_controls_EAS", "af_controls_EUR", "af_controls_MID", "beta_AFR", "beta_AMR", "beta_CSA", "beta_EAS", "beta_MID", "se_AFR", "se_AMR", "se_CSA", "se_EAS", "se_MID", "pval_AFR", "pval_AMR", "pval_CSA", "pval_EAS", "pval_MID", "low_confidence_AFR", "low_confidence_AMR", "low_confidence_CSA", "low_confidence_EAS", "low_confidence_EUR", "low_confidence_MID", "nearest_genes"]
    threads: 8
    resources:
        mem_mb = get_mem_mb,
    group: "gwas"
    script: "../scripts/misc_preprocessing.R"

# TODO rewrite pipeline to handle temporary dir, work without cd etc.
# TODO pipeline is currently handling rm of a lot of stuff
# TODO can only be run serially
rule process_gwas:
    input:
        "results/processed_gwas/{trait}_pre_pipeline.tsv.gz"
    output:
        temp_input_cp = temp("workflow/scripts/GWAS_tools/{trait}.tsv.gz"),
        processed_file = "results/processed_gwas/{trait,[^\_]+}_post_pipeline.tsv.gz"
    params:
        gwas_tools_dir = "workflow/scripts/GWAS_tools",
        pipeline_output_file = lambda w: f"workflow/scripts/GWAS_tools/{w.trait}-hg38.tsv.gz",
        temp_input_cp_name = "{trait}.tsv.gz"
    resources:
        runtime = 90
    group: "gwas"
    shell:
        """
        cp {input} {output.temp_input_cp}
        cd {params.gwas_tools_dir}

        ./pipeline_v5.3.2_beta.sh -f {wildcards.trait}.tsv.gz -b {wildcards.trait}

        cd ../../..

        mv {params.pipeline_output_file} {output.processed_file}
        """

rule recalculate_p_values:
    input:
        ancient("results/processed_gwas/{trait}_post_pipeline.tsv.gz")
    output:
        "results/processed_gwas/{trait,[^\_]+}_recalculated_p.tsv.gz"
    params:
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 8
    resources:
        mem_mb = get_mem_mb
    group: "gwas"
    script: "../scripts/recalculate_p_values.R"

rule join_pair_gwas:
    input:
        A = ancient("results/processed_gwas/{trait_A}_recalculated_p.tsv.gz"),
        B = ancient("results/processed_gwas/{trait_B}_recalculated_p.tsv.gz")
    output:
        AB = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/merged.tsv.gz"),
    threads: 8
    resources:
        mem_mb = get_mem_mb,
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
        join = 'inner',
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P',
        beta_col = 'BETA',
        se_col = 'SE',
        id_col = 'SNPID'
    group: "gwas"
    script:
        "../scripts/join_pair_gwas_stats.R"

rule make_plink_ranges:
    input:
        ("resources/1000g/hg38/eur/qc/{variant_set}/{variant_type}/chr%d.bim" % x for x in range(1,23)),
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/merged.tsv.gz"
    output:
        ("results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/matching_ids/chr%d.txt" % x for x in range(1,23))
    params:
        input_dir = "resources/1000g/hg38/eur/qc/all",
        output_dir = "results/merged_gwas/{trait_A}_and_{trait_B}/{variant_set}/{variant_type}/matching_ids",
        bim_regex = "chr%d.bim",
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT'
    threads: 8
    resources:
        mem_mb = get_mem_mb
    group: "gwas"
    script: '../scripts/make_plink_ranges.R'
