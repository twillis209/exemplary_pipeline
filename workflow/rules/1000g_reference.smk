import pandas as pd
import re

rule download_1000g_hg38_manifest:
    output:
        temp("resources/1000g/hg38/manifest.tsv")
    group: "1000g"
    shell:
        "wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv"

rule process_1000g_hg38_manifest:
    input:
        "resources/1000g/hg38/manifest.tsv"
    output:
        "resources/1000g/hg38/processed_manifest.tsv"
    group: "1000g"
    run:
        daf = pd.read_csv(input[0], sep = '\t', names = ['File', 'Byte', 'Checksum'])

        daf = daf[daf.File.str.match('.+\.vcf\.gz$')]

        daf = daf.assign(Chr=daf.File.str.extract('.+chr(\w+)\.filtered.+'))

        daf.to_csv(output[0], sep = '\t', index = False)

rule download_1000g_hg38_genotype_data:
    input:
        "resources/1000g/hg38/processed_manifest.tsv"
    output:
        temp("resources/1000g/hg38/{chr}.vcf.gz")
    resources:
        runtime = 60
    group: "1000g"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O resources/1000g/hg38/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
        else:
            shell("wget -O resources/1000g/hg38/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")

        manifest = pd.read_csv(input[0], sep = '\t')

        md5sum_expected = manifest[manifest.Chr == wildcards.chr.replace('chr', '')].Checksum.values[0]

        #shell("md5sum {output}")
        md5sum_actual = shell("md5sum {output}", read = True).split(' ')[0]

        if md5sum_expected != md5sum_actual:
            raise Exception("md5sums do not match")

rule download_1000g_hg38_sample_metadata:
     output:
        "resources/1000g/hg38/ped.txt"
     group: "1000g"
     shell:
         "wget -O resources/1000g/hg38/ped.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule vcf_to_bed:
    input:
        "resources/1000g/{assembly}/{chr}.vcf.gz"
    output:
        "resources/1000g/{assembly}/{chr}.bed",
        "resources/1000g/{assembly}/{chr}.bim",
        "resources/1000g/{assembly}/{chr}.fam"
    params:
        out = "resources/1000g/{assembly}/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1000g"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --vcf {input} --make-bed --out {params.out} --set-all-var-ids @:#:\$r:\$a --max-alleles 2 --new-id-max-allele-len 20 'truncate'"

rule make_hg38_fam_files:
     input:
         "resources/1000g/hg38/ped.txt"
     output:
         eur = "resources/1000g/hg38/eur.fam",
         afr = "resources/1000g/hg38/afr.fam",
         amr = "resources/1000g/hg38/amr.fam",
         eas = "resources/1000g/hg38/eas.fam",
         sas = "resources/1000g/hg38/sas.fam"
     script:
        "../scripts/get_fam_files.R"

rule get_ancestry_specific_samples:
     input:
        "resources/1000g/{assembly}/{chr}.bed",
        "resources/1000g/{assembly}/{chr}.bim",
        "resources/1000g/{assembly}/{chr}.fam",
        fam_file = "resources/1000g/{assembly}/{ancestry}.fam"
     output:
        temp("resources/1000g/{assembly}/{ancestry}/{chr}.bed"),
        temp("resources/1000g/{assembly}/{ancestry}/{chr}.bim"),
        temp("resources/1000g/{assembly}/{ancestry}/{chr}.fam")
     params:
        bfile = "resources/1000g/{assembly}/{chr}",
        out = "resources/1000g/{assembly}/{ancestry}/{chr}"
     threads: 8
     resources:
        mem_mb=get_mem_mb
     group: "1000g"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --keep {input.fam_file} --make-bed --silent --out {params.out}"

rule qc:
     input:
         ancient("resources/1000g/{assembly}/{ancestry}/{chr}.bed"),
         ancient("resources/1000g/{assembly}/{ancestry}/{chr}.bim"),
         ancient("resources/1000g/{assembly}/{ancestry}/{chr}.fam")
     output:
         "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}.bed",
         "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}.bim",
         "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}.fam"
     params:
        bfile = "resources/1000g/{assembly}/{ancestry}/{chr}",
        out = "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}",
     threads: 8
     resources:
        mem_mb=get_mem_mb
     group: "1000g"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --geno 0.1 --mind 0.1 --maf 0.005 --hwe 1e-50 --rm-dup 'force-first' --make-bed --silent --out {params.out}"

rule make_mhc_range:
    output:
        temp("resources/1000g/{assembly}/{ancestry}/qc/sans_mhc/mhc_range.txt")
    shell:
        """
        echo -e "6\t24000000\t45000000\tMHC" >>{output}
        """

rule remove_mhc:
    input:
        bed = "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}.bed",
        bim = "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}.bim",
        fam = "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}.fam",
        range = "resources/1000g/{assembly}/{ancestry}/qc/sans_mhc/mhc_range.txt"
    output:
        bed = "resources/1000g/{assembly}/{ancestry}/qc/sans_mhc/{chr}.bed",
        bim = "resources/1000g/{assembly}/{ancestry}/qc/sans_mhc/{chr}.bim",
        fam = "resources/1000g/{assembly}/{ancestry}/qc/sans_mhc/{chr}.fam",
    params:
        bfile = "resources/1000g/{assembly}/{ancestry}/qc/all/{chr}",
        out = "resources/1000g/{assembly}/{ancestry}/qc/sans_mhc/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1000g"
    shell:
        """
        if [ {wildcards.chr} == "chr6" ]; then
            plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --exclude 'range' {input.range} --make-bed --silent --out {params.out}
        else
            cp {input.bed} {output.bed};
            cp {input.bim} {output.bim};
            cp {input.fam} {output.fam};
        fi
        """

rule exclude_non_snps:
    input:
        bed = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{chr}.bed",
        bim = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{chr}.bim",
        fam = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{chr}.fam",
    output:
        bed = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/snps_only/{chr}.bed",
        bim = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/snps_only/{chr}.bim",
        fam = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/snps_only/{chr}.fam",
    params:
        bfile = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{chr}",
        out = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/snps_only/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1000g"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.bfile} --snps-only 'just-acgt' --make-bed --silent --out {params.out}"

rule create_pruned_ranges:
    input:
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}.bed",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}.bim",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}.fam"
    output:
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}.prune.in",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}.prune.out"
    params:
        in_bfile_stem = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}",
        r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
        prune_out = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1000g"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_bfile_stem} --indep-pairwise {wildcards.window_size} 1 {params.r2} --out {params.prune_out}"

rule prune_snps:
    input:
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}.bed",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}.bim",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}.fam",
        range_file = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}.prune.out"
    output:
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}.bed",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}.bim",
        "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}.fam"
    params:
        in_bfile_stem = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/{chr}",
        out_bfile_stem = "resources/1000g/{assembly}/{ancestry}/qc/{variant_set}/{variant_type}/pruned/{window_size}_1_{r2}/{chr}"
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1000g"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_bfile_stem} --exclude {input.range_file} --make-bed --out {params.out_bfile_stem}"
