library(data.table)

if(grepl('hg38', snakemake@input[[1]])) {
  dat <- fread(snakemake@input[[1]], header = T)

  for(x in c('eur', 'afr', 'amr', 'eas', 'sas')) {
    fwrite(dat[Superpopulation == toupper(x) & FatherID == '0' & MotherID == '0', .(0, FamilyID, SampleID, FatherID, MotherID, Sex, Phenotype = -9)], file = snakemake@output[[x]], sep = '\t', col.names = F, row.names = F, quote = F)
  }
} else if(grepl('hg19', snakemake@input[[1]])) {
  ped_dat <- fread(snakemake@input[['ped']], header = T, sep = '\t')
  
  panel_dat <- fread(snakemake@input[['panel']], header = F)

  names(panel_dat) <- c('ID', 'Population', 'Superpopulation', 'Sex')

  merged_dat <- merge(panel_dat, ped_dat[, c('Individual ID', 'Paternal ID', 'Maternal ID')], by.x = 'ID', by.y = 'Individual ID', all.x = T)

  for(x in c('eur', 'afr', 'amr', 'eas', 'sas')) {
    fwrite(merged_dat[Superpopulation == toupper(x) & `Paternal ID` == '0' & `Maternal ID` == '0', .(FamilyID = 0, ID, `Paternal ID`, `Maternal ID`, Sex, Phenotype = -9)], file = snakemake@output[[x]], sep = ' ', col.names = F, row.names = F, quote = F)
  }

  fwrite(merged_dat[`Paternal ID` == '0' & `Maternal ID` == '0', .(FamilyID = 0, ID, `Paternal ID`, `Maternal ID`, Sex, Phenotype = -9)], file = snakemake@output[['all']], sep = ' ', col.names = F, row.names = F, quote = F)
} else {
  stop("Input file contains neither hg38 nor hg19.")
}
