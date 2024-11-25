import json
import os


full_dict = {
    "input": [
        snakemake.input.r1,
        snakemake.input.r2,
    ],
    "output": [
        snakemake.params.undetermined_r1,
        snakemake.params.undetermined_r2,
        snakemake.params.undetermined_bc,
    ],
    "cellular": [
        {
            "algorithm": "mdd",
            "base": "round_A",
            "codec": {},
            "distance tolerance": [snakemake.params.dist_tolerance],
            "transform": {
                "token": snakemake.params.split_position_token,
                "knit": [snakemake.params.knit_token],
            },
        },
    ],
    # 1:0:45 - T is also included in barcode
    # 1 is segment position int, 0 is exclusive, 45 is inclusive (opposite to pheniqs documentation???)
    # maybe open issue on pheniqs
    "template": {"transform": {"token": ["0::", "1:45:", "1:0:45"]}},
}

# care for cell numbers larger than 99999 - only 5 decimal places
with open(snakemake.input.filtered_bc_list, "rt") as A:
    ls_A = list(line.strip() for line in A)
for idx, seq in enumerate(ls_A):
    full_dict["cellular"][0]["codec"][f"a{idx:05}"] = {
        "barcode": [f"{seq}"],
        "concentration": 1,
        "output": [
            f"{snakemake.params.outdir}/{idx:05}_r1.fastq.gz",
            f"{snakemake.params.outdir}/{idx:05}_r2.fastq.gz",
            f"{snakemake.params.outdir_bc}/{idx:05}_bc.fastq.gz",
        ],
    }

with open(snakemake.output.split_pheniqs_config, "w") as fp:
    json.dump(full_dict, fp)
