wildcard_constraints:
    input_name = "[A-Za-z0-9\+-]+"

extraneous_columns = ["af_cases_meta_hq" , "af_controls_meta_hq", "beta_meta_hq", "se_meta_hq", "pval_meta_hq", "pval_heterogeneity_hq", "af_cases_meta", "af_controls_meta", "beta_meta", "se_meta", "pval_meta", "pval_heterogeneity", "af_cases_AFR", "af_cases_AMR", "af_cases_CSA", "af_cases_EAS", "af_cases_EUR", "af_cases_MID", "af_controls_AFR", "af_controls_AMR", "af_controls_CSA", "af_controls_EAS", "af_controls_EUR", "af_controls_MID", "beta_AFR", "beta_AMR", "beta_CSA", "beta_EAS", "beta_MID", "se_AFR", "se_AMR", "se_CSA", "se_EAS", "se_MID", "pval_AFR", "pval_AMR", "pval_CSA", "pval_EAS", "pval_MID", "low_confidence_AFR", "low_confidence_AMR", "low_confidence_CSA", "low_confidence_EAS", "low_confidence_EUR", "low_confidence_MID", "nearest_genes", "neglog10_pval_meta_hq", "neglog10_pval_heterogeneity_hq", "neglog10_pval_meta", "neglog10_pval_heterogeneity", "neglog10_pval_AFR", "neglog10_pval_CSA", "neglog10_pval_EAS"]

rule rehead:
    input:
        "resources/gwas/{input_name}.tsv.gz"
    output:
        temp("results/gwas_pipeline/reheaded/{input_name}.tsv.gz")
    params:
        columns_to_drop = extraneous_columns,
        pan_ukb_p_column = 'pval_EUR',
        pan_ukb_neglog10_p_column = 'neglog10_pval_EUR',
        pan_ukb_beta_column = 'beta_EUR',
        pan_ukb_se_column = 'se_EUR'
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/rehead.R"


rule fix_alleles_and_id:
    input:
        "results/gwas_pipeline/reheaded/{input_name}.tsv.gz"
    output:
        temp("results/gwas_pipeline/reheaded/fixed_alleles/{input_name}.tsv.gz")
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/fix_alleles_and_id.R"

rule check_for_minimal_column_set:
    input:
        "results/gwas_pipeline/reheaded/fixed_alleles/{input_name}.tsv.gz"
    output:
        temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/{input_name}.tsv.gz")
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/check_for_minimal_column_set.R"

rule recalculate_missing_summary_statistics:
    input:
        "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/{input_name}.tsv.gz"
    output:
        temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/{input_name}.tsv.gz")
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/recalculate_missing_sumstats.R"

rule detect_build:
    input:
        manifest = "resources/gwas_pipeline/build_manifest.tsv",
        sumstats = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/{input_name}.tsv.gz"
    output:
        temp("results/gwas_pipeline/build_detection/{input_name}.tsv")
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/detect_build.R"

rule create_bed_file_for_liftover:
    input:
        "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/{input_name}.tsv.gz"
    output:
        temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}.bed")
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/create_bedfile.R"

rule prepare_file_for_liftover:
    input:
        sumstats = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/{input_name}.tsv.gz",
        build = "results/gwas_pipeline/build_detection/{input_name}.tsv"
    output:
        prepared = temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_prepared.tsv.gz"),
        build = temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_build.txt")
    threads: 8
    group: "gwas"
    script: "../scripts/gwas_pipeline/prepare_file_for_liftover.R"

rule liftover:
    input:
        build = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_build.txt",
        bed_file = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}.bed",
        hg19_chainfile = "resources/gwas_pipeline/hg19ToHg38.over.chain.gz"
    output:
        lifted = temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_lifted.bed"),
        unlifted = temp("results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_unlifted.bed")
    threads: 1
    group: "gwas"
    run:
        with open(input.build, 'r') as fh:
            assembly = fh.readline().strip()

        if assembly == 'hg19':
            shell("liftOver {input.bed_file} {input.hg19_chainfile} {output.lifted} {output.unlifted}")
        elif assembly == 'hg38':
            print("Already in hg38, copying...")
            shell("cp {input.bed_file} {output.lifted}")
            shell("touch {output.unlifted}")
        else:
            raise Exception('Can only liftover from hg19 at the moment')

rule merge_liftovered_rows_with_summary_statistics:
    input:
        lifted = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_lifted.bed",
        sumstats = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_prepared.tsv.gz",
        build = "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/{input_name}_build.txt"
    threads: 8
    output:
        "results/processed_gwas/{input_name}.tsv.gz"
    group: "gwas"
    script: "../scripts/gwas_pipeline/merge_liftovered_rows_with_sumstats.R"
