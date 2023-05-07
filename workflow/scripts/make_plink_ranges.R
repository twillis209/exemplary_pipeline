library(data.table)
options(warn = 2)

setDTthreads(snakemake@threads)

gwas_file <- snakemake@input[['gwas_file']]

input_dir <- snakemake@params[['input_dir']]
output_dir <- snakemake@params[['output_dir']]
bim_regex <- snakemake@params[['bim_regex']]
chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
tmpdir <- snakemake@resources[['tmpdir']]

gwas_dat <- fread(gwas_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col), tmpdir = tmpdir)

gwas_dat[ , chr_col := as.character(get(chr_col))]

for(i in 1:22) {
  bim_dat <- fread(file.path(input_dir, sprintf(bim_regex, i)), sep = '\t', header = F, col.names = c('CHR38', 'ID', 'Cm', 'BP38', 'A1', 'A2'))

  bim_dat[, CHR38 := as.character(CHR38)]

  bim_join <- merge(bim_dat, gwas_dat, by.x = c(chr_col, bp_col), by.y = c('CHR38', 'BP38'), sort = F)

  # Make sure alleles match
  bim_join <- bim_join[(get(ref_col) == A1 & get(alt_col) == A2) | (get(ref_col) == A2 & get(alt_col) == A1)]

  bim_join <- bim_join[, .(ID)]

  fwrite(bim_join, file = file.path(output_dir, sprintf("chr%d.txt", i)), row.names = F, sep = ' ', col.names = F, quote = F)
}
