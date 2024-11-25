import pandas as pd

user_input = int(snakemake.config["read_threshold"])

df = pd.read_csv(snakemake.input.observed_bc_list, sep="\s+", names=["number", "bc"])
df1 = df[(df["number"] >= user_input) & (~df["bc"].str.contains("="))].sort_values(
    "number", ascending=False
)

df1.head(snakemake.config["cell_threshold"]).bc.str.replace(
    "-", ""
).to_csv(snakemake.output.filtered_bc_list, sep="\n", header=False, index=False)
