import glob
import os
# import subprocess
# from pathlib import Path
# import numpy as np
import pandas as pd
# from io import StringIO


coverage_file_list = snakemake.input.coverage_tables

# filter out undetermined
filtered_list = [item for item in coverage_file_list  if "undetermined" not in item]


df_list = []

for table_pth in filtered_list:


    bc_num = table_pth[-18:-13] # minimap TODO:remove naming altogether

    df = pd.read_csv(table_pth, sep='\t')
    df['barcode_num'] = f"{bc_num}"

    df_list.append(df)

result_df = pd.concat(df_list)

result_df.to_csv(snakemake.output.single_coverage_table)
