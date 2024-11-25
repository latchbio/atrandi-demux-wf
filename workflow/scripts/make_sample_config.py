import json
from pathlib import Path
import os


full_dict = {
    "input": [snakemake.input.sampled_fastq],
    "output": [snakemake.params.demux_sampled_bam],
    "cellular": [
        {
            "algorithm": "mdd",
            "base": "round_D",
            "codec": {},
            "distance tolerance": [snakemake.params.dist_tolerance],
            "transform": {"token": [snakemake.params.token_D]},
        },
        {
            "algorithm": "mdd",
            "base": "round_C",
            "codec": {},
            "distance tolerance": [snakemake.params.dist_tolerance],
            "transform": {"token": [snakemake.params.token_C]},
        },
        {
            "algorithm": "mdd",
            "base": "round_B",
            "codec": {},
            "distance tolerance": [snakemake.params.dist_tolerance],
            "transform": {"token": [snakemake.params.token_B]},
        },
        {
            "algorithm": "mdd",
            "base": "round_A",
            "codec": {},
            "distance tolerance": [snakemake.params.dist_tolerance],
            "transform": {"token": [snakemake.params.token_A]},
        },
    ],
    "template": {"transform": {"token": ["0::"]}},
}


# care not to mix 0,1,2 and 3 position values
with open(snakemake.input.bcD, "rt") as D:
    ls_D = list(line.strip() for line in D)
for idx, seq in enumerate(ls_D):
    full_dict["cellular"][0]["codec"][f"d{idx:02}"] = {
        "barcode": [f"{seq}"],
        "concentration": 1,
    }
with open(snakemake.input.bcC, "rt") as C:
    ls_C = list(line.strip() for line in C)
for idx, seq in enumerate(ls_C):
    full_dict["cellular"][1]["codec"][f"c{idx:02}"] = {
        "barcode": [f"{seq}"],
        "concentration": 1,
    }
with open(snakemake.input.bcB, "rt") as B:
    ls_B = list(line.strip() for line in B)
for idx, seq in enumerate(ls_B):
    full_dict["cellular"][2]["codec"][f"b{idx:02}"] = {
        "barcode": [f"{seq}"],
        "concentration": 1,
    }
with open(snakemake.input.bcA, "rt") as A:
    ls_A = list(line.strip() for line in A)
for idx, seq in enumerate(ls_A):
    full_dict["cellular"][3]["codec"][f"a{idx:02}"] = {
        "barcode": [f"{seq}"],
        "concentration": 1,
    }

with open(snakemake.output.sample_pheniqs_config, "w") as fp:
    json.dump(full_dict, fp)
