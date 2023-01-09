import numpy as np
import pandas as pd
import logging
from easyfinemap.loci import Loci
from easyfinemap.constant import ColName
from easyfinemap.utils import get_significant_snps, make_SNPID_unique
from easyfinemap import LDRef
from easyfinemap.sumstat import SumStat
from easyfinemap.easyfinemap import EasyFinemap

from pathlib import Path

# l = Loci()
# df = pd.read_csv('./exampledata/noEAF_noMAF.txt.gz', sep='\t')
# leadsnp, indep_loci = l.identify_indep_loci(
#     df,
#     method='conditional',
#     ldref='./exampledata/LDREF/EUR.valid.chr{chrom}',
#     sample_size=1000,
#     use_ref_EAF=True,
#     threads=2,
# )

df = pd.read_csv('./exampledata/noEAF_noMAF.txt.gz',sep='\t')
df = SumStat(df)
df = df.standarize()

indep_loci = pd.read_csv('exampledata/conditional.loci.txt', sep='\t')

chrom, start, end = indep_loci.loc[0, ['CHR', 'START', 'END']]
sub_df = df.loc[(df[ColName.CHR] == chrom) & (df[ColName.BP] >= start) & (df[ColName.BP] <= end)].copy()
ld = LDRef()
sub_df_ol = ld.intersect(sub_df, f'./exampledata/LDREF/EUR.valid.chr{chrom}', './tmp/intersc', use_ref_EAF=True)
ld.make_ld('./tmp/intersc', './tmp/ld')
ef = EasyFinemap()
finemap_pp = ef.run_paintor(sub_df_ol, ld_matrix='./tmp/ld.ld')
