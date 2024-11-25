from dataclasses import dataclass
from typing import List

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata.latch import LatchAuthor
from latch.types.metadata.snakemake import SnakemakeParameter
from latch.types.metadata.snakemake_v2 import SnakemakeV2Metadata


@dataclass
class Sample:
    name: str
    r1: LatchFile
    r2: LatchFile


metadata = SnakemakeV2Metadata(
    display_name="Atrandi Pheniqs Demultiplexing",
    author=LatchAuthor(),
    parameters={
        "samples": SnakemakeParameter(
            display_name="Samples",
            type=List[Sample],
            section_title="Folders",
            samplesheet=True,
        ),
        "results_dir": SnakemakeParameter(
            display_name="Results Dir",
            type=LatchOutputDir,
            default=LatchDir("latch://1721.account/atrandi-results"),
        ),
        "bc_D": SnakemakeParameter(
            display_name="Bc D",
            type=LatchFile,
            default=LatchFile("latch://1721.account/Barcodes/bcD_24.txt"),
            hidden=True,
            section_title="Barcodes",
        ),
        "bc_C": SnakemakeParameter(
            display_name="Bc C",
            type=LatchFile,
            default=LatchFile("latch://1721.account/Barcodes/bcC_24.txt"),
            hidden=True,
        ),
        "bc_B": SnakemakeParameter(
            display_name="Bc B",
            type=LatchFile,
            default=LatchFile("latch://1721.account/Barcodes/bcB_24.txt"),
            hidden=True,
        ),
        "bc_A": SnakemakeParameter(
            display_name="Bc A",
            type=LatchFile,
            default=LatchFile("latch://1721.account/Barcodes/bcA_24.txt"),
            hidden=True,
        ),
        "random_seed": SnakemakeParameter(
            display_name="Random Seed",
            type=int,
            default=42,
            section_title="General",
        ),
        "fraction_to_sample": SnakemakeParameter(
            display_name="Read number to sample",
            type=int,
            default=1000000,
        ),
        "read_threshold": SnakemakeParameter(
            display_name="Read Threshold",
            type=int,
            default=1,
        ),
        "cell_threshold": SnakemakeParameter(
            display_name="Cell Threshold",
            type=int,
            default=100000,
        ),
        "threads": SnakemakeParameter(
            display_name="Threads",
            type=int,
            default=16,
        ),
        "sample_dist_tolerance": SnakemakeParameter(
            display_name="Sample Dist Tolerance",
            type=int,
            default=1,
        ),
        "split_dist_tolerance": SnakemakeParameter(
            display_name="Split Dist Tolerance",
            type=int,
            default=1,
        ),
        "sample_token_D": SnakemakeParameter(
            display_name="Sample Token D",
            type=str,
            default="0:0:8",
            hidden=True,
            section_title="Tokens",
        ),
        "sample_token_C": SnakemakeParameter(
            display_name="Sample Token C",
            type=str,
            default="0:12:20",
            hidden=True,
        ),
        "sample_token_B": SnakemakeParameter(
            display_name="Sample Token B",
            type=str,
            default="0:24:32",
            hidden=True,
        ),
        "sample_token_A": SnakemakeParameter(
            display_name="Sample Token A",
            type=str,
            default="0:36:44",
            hidden=True,
        ),
        "split_position_token": SnakemakeParameter(
            display_name="Split Position Token",
            type=List[str],
            default=["1:0:8", "1:12:20", "1:24:32", "1:36:44"],
            hidden=True,
        ),
        "split_knit_token": SnakemakeParameter(
            display_name="Split Knit Token",
            type=str,
            default="0:1:2:3",
            hidden=True,
        ),
    },
)
