wildcard_constraints:
    snp_set = 'with_mhc|sans_mhc',
    chr = "chr[0-9XY]{1,2}",
    assembly = "hg19|hg38",
    ancestry = "eur|afr|amr|eas|sas|all",
    variant_set = "all|sans_mhc",
    variant_type = "all|snps_only",
    window_size = "\\d+kb",
    r2 = "0_\\d+"
