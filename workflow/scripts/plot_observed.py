import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# todo - clean up/organise this

df = pd.read_csv(
    snakemake.input.observed_bc_list, sep="\s+", names=["read_count_per_bc", "bc_name"]
)
df1 = df[~df["bc_name"].str.contains("=")].reset_index(drop=True)

barcodes_reads = df1.sort_values("read_count_per_bc", ascending=False)

sns.set_style("whitegrid")


bins = 25
figw = 12
figh = 10

data = barcodes_reads.read_count_per_bc

xmin = min(data)
xmax = max(data) * 10.0
ts = snakemake.config["read_threshold"]

# since log, ts = 1 can't be plotted
if ts <= 1:
    ts_plot = 3
else:
    ts_plot = ts

print(ts_plot)

sample_name = snakemake.wildcards


def startfig(
    w=4, h=2, rows=1, columns=1, wrs=None, hrs=None, frameon=True, return_first_ax=True
):
    """
    for initiating figures, w and h in centimeters
    example of use:
    a,fig,gs = startfig(w=10,h=2.2,rows=1,columns=3,wr=[4,50,1],hrs=None,frameon=True)
    hrs - height ratios
    wrs - width ratios
    frameon - whether first axes with frame

    returns:
    if return_first_ax=True
    a,fig,gs
    else
    fig,gs
    """

    ratio = 0.393701  # 1 cm in inch
    myfigsize = (w * ratio, h * ratio)
    fig = plt.figure(figsize=(myfigsize))
    gs = mpl.gridspec.GridSpec(rows, columns, width_ratios=wrs, height_ratios=hrs)
    if return_first_ax == True:
        a = fig.add_subplot(gs[0, 0], frameon=frameon)
        return a, fig, gs
    else:
        return fig, gs


a, fig, gs = startfig(figw, figh, rows=1)

a.hist(
    np.log10(data),
    bins=bins,
    lw=1,
    facecolor="0.5",
    edgecolor="k",
    histtype="stepfilled",
)
a.set_xlabel("log10(# reads)")
a.set_ylabel("# barcodes")
# plotting ts_plot instead of ts
a.axvline(np.log10(ts_plot), color="r")
a.set_title("%s\nBarcode distribution" % sample_name)

m = data >= ts  # mask where data is more or equal to ts
abs_above = m.sum()
frac_above = m.sum() / len(m)
abs_below = (~m).sum()
frac_below = 1.0 - frac_above

a.text(
    0,
    -0.3,
    "threshold = %d\n%d (%.1f%%) barcodes >= threshold\n%d (%.1f%%) barcodes < threshold"
    % (ts, abs_above, frac_above * 100.0, abs_below, frac_below * 100.0),
    horizontalalignment="left",
    verticalalignment="top",
    transform=plt.gca().transAxes,
    bbox=dict(facecolor="1", alpha=0.8),
)


a.set_xlim(np.log10(xmin), np.log10(xmax))
gs.tight_layout(fig)
plt.savefig(snakemake.output.observed_barcode_hist)


a, fig, gs = startfig(figw, figh, rows=1)

a.hist(
    np.log10(data),
    bins=bins,
    weights=data,
    lw=1,
    facecolor="c",
    edgecolor="k",
    histtype="stepfilled",
)
a.set_xlabel("log10(# reads)")
a.set_ylabel("# reads coming from bin")
# plotting ts_plot instead of ts
a.axvline(np.log10(ts_plot), color="r", label="threshold = %d" % ts)
a.set_title("%s\nRead distribution" % sample_name)

m = data >= ts
abs_above = data[m].sum()
frac_above = abs_above / data.sum()
abs_below = data[~m].sum()
frac_below = 1.0 - frac_above

a.text(
    0,
    -0.3,
    "threshold = %d\n%d (%.1f%%) reads >= threshold\n%d (%.1f%%) reads < threshold"
    % (ts, abs_above, frac_above * 100.0, abs_below, frac_below * 100.0),
    horizontalalignment="left",
    verticalalignment="top",
    transform=plt.gca().transAxes,
    bbox=dict(facecolor="1", alpha=0.8),
)

a.set_xlim(np.log10(xmin), np.log10(xmax))
gs.tight_layout(fig)
plt.savefig(snakemake.output.observed_read_hist)
