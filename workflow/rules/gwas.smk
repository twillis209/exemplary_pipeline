def get_url(w):
    url = metadata_daf.loc[metadata_daf['abbrv'] == w.trait, 'URL'].values[0]

    return url

def get_expected_md5sum(w):
    md5sum = metadata_daf.loc[metadata_daf['abbrv'] == w.trait, 'md5'].values[0]

    return md5sum

rule download_gwas:
    output:
        "resources/gwas/{trait}.tsv.gz"
    params:
        url = get_url,
        expected_md5sum = get_expected_md5sum
    # NB: We parallelise compression where necessary
    threads: lambda w: 8 if 'tsv.bgz' in get_url(w) else 1
    resources:
        runtime = 15
    group: "gwas"
    run:
        if '.tsv.bgz' in params.url:
            temp_bgz_output = output[0].replace('tsv.gz', 'tsv.bgz')
            shell(f"wget -O {temp_bgz_output} {params.url}")

            actual_md5sum = shell(f"md5sum {temp_bgz_output}", read = True).split(' ')[0]

            print(params.expected_md5sum)
            print(actual_md5sum)

            if params.expected_md5sum != actual_md5sum:
                raise Exception(f"md5sums do not match for {wildcards.trait}")
            else:
                print("md5sum checked")

            temp_tsv_output = output[0].replace('tsv.gz', 'tsv')
            shell(f"bgzip -@ {threads} -d {temp_bgz_output}")
            shell(f"pigz -p {threads} {temp_tsv_output}")
        # TODO elif handle T1D, which is apparently not compressed
        else:
            shell("wget -O {output} {params.url}")

            actual_md5sum = shell("md5sum {output}", read = True).split(' ')[0]

            print(params.expected_md5sum)
            print(actual_md5sum)

            if params.expected_md5sum != actual_md5sum:
                raise Exception(f"md5sums do not match for {wildcards.trait}")
            else:
                print("md5sum checked")

rule join_pair_gwas:
    input:
        A = "results/processed_gwas/{trait_A}.tsv.gz",
        B = "results/processed_gwas/{trait_B}.tsv.gz"
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
    conda: "../envs/exemplary_pipeline.yaml"
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
    conda: "../envs/exemplary_pipeline.yaml"
    script: '../scripts/make_plink_ranges.R'
