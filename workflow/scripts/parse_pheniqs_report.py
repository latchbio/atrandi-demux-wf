import pandas as pd
import json

with open(snakemake.input.split_report, "r") as file:
    json_data = json.load(file)


def extract_barcode_info(data):
    barcode_info = []
    for entry in data.get("cellular", []):
        for classified_entry in entry.get("classified", []):
            barcode_info.append(
                {
                    "barcode": classified_entry["barcode"][0],
                    "count": classified_entry["count"],
                    "fastq_index": classified_entry["index"],
                }
            )
    return barcode_info


barcode_data = extract_barcode_info(json_data)
df = pd.DataFrame(barcode_data)
df["fastq_index"] = df["fastq_index"].apply(lambda x: str(x - 1).zfill(5))
df.to_csv(snakemake.output.split_barcode_table, header=False, index=False)
