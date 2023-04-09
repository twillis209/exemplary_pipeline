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

rule join_pair_gwas:
    input:
        A = ancient("results/processed_gwas/{trait_A}.tsv.gz"),
        B = ancient("results/processed_gwas/{trait_B}.tsv.gz")
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
