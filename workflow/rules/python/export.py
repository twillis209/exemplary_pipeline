import os
import pandas as pd
from scipy.stats import chi2

def compile_sumher_files(input_files, output_file):
    d = []

    for x in input_files:
        trait_A, trait_B = x.split('/')[3].split('_')
        snp_set = x.split('/')[4]

        with open(x, 'r') as infile:
            line = infile.readline()
            line = infile.readline()

        # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
        _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()

        rg_z = float(rg)/float(rg_se)

        rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)

        d.append(
            {
                'trait.A' : trait_A,
                'trait.B' : trait_B,
                'snp.set' : snp_set,
                'h2.A.obs.sr' : float(h2_A),
                'h2.A.obs.se.sr' : float(h2_A_se),
                'h2.B.obs.sr' : float(h2_B),
                'h2.B.obs.se.sr' : float(h2_B_se),
                'gcov.obs.sr' : float(cov),
                'gcov.obs.se.sr' : float(cov_se),
                'rg.sr' : float(rg),
                'rg.se.sr' : float(rg_se),
                'rg.z.sr' : rg_z,
                'rg.p.sr' : rg_p
            }
        )

        print(x)

    pd.DataFrame(d).to_csv(output_file, sep = '\t', index = False)

if __name__ == '__main__':
    compile_sumher_files(snakemake.input, snakemake.output[0])
