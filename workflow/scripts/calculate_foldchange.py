import pandas as pd
import numpy as np
import scipy
import scipy.stats


df = pd.read_csv(snakemake.input.coverage_table,index_col=0, dtype = {'barcode_num': str})
df_lookup = pd.read_csv(snakemake.input.lookup_table,index_col=None, header=None, names = ["bc","read_count","file_num"], dtype = {'file_num': str})

lookup_dict = dict(zip(df_lookup['file_num'], df_lookup['bc']))
lookup_dict_2 = dict(zip(df_lookup['file_num'], df_lookup['read_count']))


df['barcode'] = df['barcode_num'].map(lookup_dict)
df['read_count_per_bc'] = df['barcode_num'].map(lookup_dict_2)

# filter some chromosomes - mixed human and bacterial
# add all possible names - check by star index!!!!

# technical debt starting to pile up - split this to separate rules for foldchange, filtering and lookup

if snakemake.params.filter_chromosomes:
    filter_list = snakemake.params.filter_chromosomes
    df = df[df["#rname"].isin(filter_list)]
    

# filter_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT","NZ_LR881938.1","NC_000964.3", "Bacillus_NZ_CP034943.1", "E_Coli_NC_000913.3" ]
# df = df[df["#rname"].isin(filter_list)]


# calculate fold change
def fold_change(row):
    maxcovs = 1.-scipy.stats.poisson.pmf(0,row["meandepth"])
    # maxcovs = maxcovs*100.
    return maxcovs

df['fold_change'] = df.apply(fold_change, axis=1)

df.to_csv(snakemake.output.fold_change_table)