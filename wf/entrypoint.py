import json
import os
import shutil
import subprocess
import sys
import typing
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from pathlib import Path
from urllib.parse import urlparse

import requests
from latch.resources.tasks import custom_task, snakemake_runtime_task
from latch.resources.workflow import workflow
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.services.register.utils import import_module_by_path

import_module_by_path(Path("latch_metadata/__init__.py"))

import latch.types.metadata.snakemake_v2 as smv2


def get_config_val(val):
    if isinstance(val, list):
        return [get_config_val(x) for x in val]
    if isinstance(val, dict):
        return {k: get_config_val(v) for k, v in val.items()}
    if isinstance(val, (LatchFile, LatchDir)):
        if val.remote_path is None:
            return str(val.path)

        parsed = urlparse(val.remote_path)
        domain = parsed.netloc
        if domain == "":
            domain = "inferred"

        return f"/ldata/{domain}{parsed.path}"

    if isinstance(val, (int, float, bool, type(None))):
        return val
    if is_dataclass(val):
        return {f.name: get_config_val(getattr(val, f.name)) for f in fields(val)}
    if isinstance(val, Enum):
        return val.value

    return str(val)


@dataclass
class Sample:
    name: str
    r1: LatchFile
    r2: LatchFile


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 24 * 30,
            "version": 2,
            "snakemake": True,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@snakemake_runtime_task(cpu=1, memory=2, storage_gib=50)
def snakemake_runtime(
    pvc_name: str,
    samples: typing.List[Sample],
    results_dir: LatchOutputDir = LatchDir("latch://36970.account/results"),
    bc_D: LatchFile = LatchFile("latch://36970.account/barcodes/bcD_24.txt"),
    bc_C: LatchFile = LatchFile("latch://36970.account/barcodes/bcC_24.txt"),
    bc_B: LatchFile = LatchFile("latch://36970.account/barcodes/bcB_24.txt"),
    bc_A: LatchFile = LatchFile("latch://36970.account/barcodes/bcA_24.txt"),
    random_seed: int = 42,
    fraction_to_sample: int = 1000000,
    read_threshold: int = 1,
    cell_threshold: int = 100000,
    threads: int = 16,
    sample_dist_tolerance: int = 1,
    split_dist_tolerance: int = 1,
    sample_token_D: str = "0:0:8",
    sample_token_C: str = "0:12:20",
    sample_token_B: str = "0:24:32",
    sample_token_A: str = "0:36:44",
    split_position_token: typing.List[str] = ["1:0:8", "1:12:20", "1:24:32", "1:36:44"],
    split_knit_token: str = "0:1:2:3",
):
    print(f"Using shared filesystem: {pvc_name}")

    shared = Path("/snakemake-workdir")
    snakefile = shared / "Snakefile.demux"

    config = {
        "samples": get_config_val(samples),
        "results_dir": get_config_val(results_dir),
        "bc_D": get_config_val(bc_D),
        "bc_C": get_config_val(bc_C),
        "bc_B": get_config_val(bc_B),
        "bc_A": get_config_val(bc_A),
        "random_seed": get_config_val(random_seed),
        "fraction_to_sample": get_config_val(fraction_to_sample),
        "read_threshold": get_config_val(read_threshold),
        "cell_threshold": get_config_val(cell_threshold),
        "threads": get_config_val(threads),
        "sample_dist_tolerance": get_config_val(sample_dist_tolerance),
        "split_dist_tolerance": get_config_val(split_dist_tolerance),
        "sample_token_D": get_config_val(sample_token_D),
        "sample_token_C": get_config_val(sample_token_C),
        "sample_token_B": get_config_val(sample_token_B),
        "sample_token_A": get_config_val(sample_token_A),
        "split_position_token": get_config_val(split_position_token),
        "split_knit_token": get_config_val(split_knit_token),
    }

    print(f"Config: {json.dumps(config, indent=2)}")
    print()

    config_path = (shared / "__latch.config.json").resolve()
    config_path.write_text(json.dumps(config, indent=2))

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        ".snakemake",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    cmd = [
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--configfile",
        str(config_path),
        "--executor",
        "latch",
        "--jobs",
        "1000",
        "--quiet",
        "--verbose",
    ]

    print("Launching Snakemake Runtime")
    print(" ".join(cmd), flush=True)

    failed = False
    try:
        subprocess.run(cmd, cwd=shared, check=True)
    except subprocess.CalledProcessError:
        failed = True
    finally:
        if not failed:
            return

        sys.exit(1)


@workflow(smv2._snakemake_v2_metadata)
def snakemake_v2_atrandi_pheniqs_demultiplexing(
    samples: typing.List[Sample],
    results_dir: LatchOutputDir = LatchDir("latch://36970.account/results"),
    bc_D: LatchFile = LatchFile("latch://36970.account/barcodes/bcD_24.txt"),
    bc_C: LatchFile = LatchFile("latch://36970.account/barcodes/bcC_24.txt"),
    bc_B: LatchFile = LatchFile("latch://36970.account/barcodes/bcB_24.txt"),
    bc_A: LatchFile = LatchFile("latch://36970.account/barcodes/bcA_24.txt"),
    random_seed: int = 42,
    fraction_to_sample: int = 1000000,
    read_threshold: int = 1,
    cell_threshold: int = 100000,
    threads: int = 16,
    sample_dist_tolerance: int = 1,
    split_dist_tolerance: int = 1,
    sample_token_D: str = "0:0:8",
    sample_token_C: str = "0:12:20",
    sample_token_B: str = "0:24:32",
    sample_token_A: str = "0:36:44",
    split_position_token: typing.List[str] = ["1:0:8", "1:12:20", "1:24:32", "1:36:44"],
    split_knit_token: str = "0:1:2:3",
):
    """
    Sample Description
    """

    snakemake_runtime(
        pvc_name=initialize(),
        samples=samples,
        results_dir=results_dir,
        bc_D=bc_D,
        bc_C=bc_C,
        bc_B=bc_B,
        bc_A=bc_A,
        random_seed=random_seed,
        fraction_to_sample=fraction_to_sample,
        read_threshold=read_threshold,
        cell_threshold=cell_threshold,
        threads=threads,
        sample_dist_tolerance=sample_dist_tolerance,
        split_dist_tolerance=split_dist_tolerance,
        sample_token_D=sample_token_D,
        sample_token_C=sample_token_C,
        sample_token_B=sample_token_B,
        sample_token_A=sample_token_A,
        split_position_token=split_position_token,
        split_knit_token=split_knit_token,
    )
