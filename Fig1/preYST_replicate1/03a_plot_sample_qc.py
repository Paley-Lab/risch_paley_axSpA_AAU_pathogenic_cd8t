"""
plot_sample_qc.py
-----------------
Reproduces the three sample QC barplots from metrics.csv.

Usage:
    python plot_sample_qc.py                        # reads metrics.csv in cwd
    python plot_sample_qc.py --input /path/to/metrics.csv --output qc_plots.pdf
    python plot_sample_qc.py --threshold 0.70       # change low-coverage threshold

Requirements:
    pip install pandas matplotlib
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from pathlib import Path


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

COLORS = {
    "mapped":    "#185FA5",
    "no_flank":  "#0F6E56",
    "mismatch":  "#BA7517",
    "no_ref":    "#A32D2D",
    "fail_bar":  "#993C1D",
    "threshold": "#E24B4A",
}


# ---------------------------------------------------------------------------
# Data loading & derived columns
# ---------------------------------------------------------------------------

def load_and_enrich(path: str, threshold: float) -> pd.DataFrame:
    df = pd.read_csv(path)

    # Support both the old format (no Total_reads) and the new enriched format
    count_cols = ["Mapped_reads", "No_flanking_sequence_hit",
                  "R1_R2_var_seqs_mismatched", "var_seq_no_ref_match"]

    # Drop any non-sample rows appended by compute_metrics (variant totals / blank rows)
    df = df.dropna(subset=["Mapped_reads"]).copy()
    df = df[df["Sample"].notna()].copy()

    if "Total_reads" not in df.columns:
        df["Total_reads"] = df[count_cols].sum(axis=1)

    df["Mapping_rate_pct"] = df["Mapped_reads"] / df["Total_reads"] * 100
    df["No_flanking_pct"]  = df["No_flanking_sequence_hit"]  / df["Total_reads"] * 100
    df["Mismatch_pct"]     = df["R1_R2_var_seqs_mismatched"] / df["Total_reads"] * 100
    df["No_ref_pct"]       = df["var_seq_no_ref_match"]       / df["Total_reads"] * 100
    df["Low_coverage"]     = df["Mapping_rate_pct"] < threshold * 100

    # Shorten long sample names for axis labels
    df["Label"] = (
        df["Sample"]
        .str.replace(r"_Rep", " R", regex=True)
    )

    return df


# ---------------------------------------------------------------------------
# Chart 1 — Mapping rate per sample
# ---------------------------------------------------------------------------

def plot_mapping_rate(ax: plt.Axes, df: pd.DataFrame, threshold_pct: float) -> None:
    bar_colors = [
        COLORS["fail_bar"] if low else COLORS["mapped"]
        for low in df["Low_coverage"]
    ]
    bars = ax.bar(df["Label"], df["Mapping_rate_pct"],
                  color=bar_colors, width=0.6, zorder=3)

    # Threshold line
    ax.axhline(threshold_pct, color=COLORS["threshold"],
               linewidth=1.5, linestyle="--", zorder=4)
    ax.text(len(df) - 0.5, threshold_pct + 1.5,
            f"{threshold_pct:.0f}% threshold",
            color=COLORS["threshold"], fontsize=9, ha="right")

    # Value labels on bars
    for bar, val in zip(bars, df["Mapping_rate_pct"]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.8,
                f"{val:.1f}%", ha="center", va="bottom", fontsize=8)

    ax.set_ylim(0, 105)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())
    ax.set_title("Mapping rate per sample", fontsize=11, fontweight="medium", pad=8)
    ax.set_ylabel("% of total reads", fontsize=9)
    _style_ax(ax, df)

    # Legend
    pass_patch = mpatches.Patch(color=COLORS["mapped"],   label="≥ threshold")
    fail_patch = mpatches.Patch(color=COLORS["fail_bar"], label="< threshold")
    ax.legend(handles=[pass_patch, fail_patch], fontsize=8, loc="lower right")


# ---------------------------------------------------------------------------
# Chart 2 — Stacked read-fate breakdown (%)
# ---------------------------------------------------------------------------

def plot_stacked_pct(ax: plt.Axes, df: pd.DataFrame) -> None:
    bottom = pd.Series([0.0] * len(df))

    layers = [
        ("Mapped",          "Mapping_rate_pct", COLORS["mapped"]),
        ("No flanking hit", "No_flanking_pct",  COLORS["no_flank"]),
        ("R1/R2 mismatch",  "Mismatch_pct",     COLORS["mismatch"]),
        ("No ref match",    "No_ref_pct",        COLORS["no_ref"]),
    ]

    for label, col, color in layers:
        ax.bar(df["Label"], df[col], bottom=bottom,
               color=color, width=0.6, label=label, zorder=3)
        bottom += df[col].values

    ax.set_ylim(0, 102)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())
    ax.set_title("Read fate breakdown per sample (% of total)", fontsize=11,
                 fontweight="medium", pad=8)
    ax.set_ylabel("% of total reads", fontsize=9)
    _style_ax(ax, df)
    ax.legend(fontsize=8, loc="lower right", ncol=2)


# ---------------------------------------------------------------------------
# Chart 3 — Absolute failure counts
# ---------------------------------------------------------------------------

def plot_absolute_failures(ax: plt.Axes, df: pd.DataFrame) -> None:
    n = len(df)
    x = range(n)
    width = 0.25
    offsets = [-width, 0, width]

    cols_colors_labels = [
        ("No_flanking_sequence_hit",  COLORS["no_flank"],  "No flanking hit"),
        ("R1_R2_var_seqs_mismatched", COLORS["mismatch"],  "R1/R2 mismatch"),
        ("var_seq_no_ref_match",      COLORS["no_ref"],    "No ref match"),
    ]

    for (col, color, label), offset in zip(cols_colors_labels, offsets):
        positions = [xi + offset for xi in x]
        ax.bar(positions, df[col], width=width,
               color=color, label=label, zorder=3)

    ax.set_xticks(list(x))
    ax.set_xticklabels(df["Label"], rotation=30, ha="right", fontsize=9)
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{v / 1e6:.1f}M")
    )
    ax.set_title("Absolute read counts per failure category", fontsize=11,
                 fontweight="medium", pad=8)
    ax.set_ylabel("Read count", fontsize=9)
    _style_ax(ax, df, set_xticks=False)
    ax.legend(fontsize=8, loc="upper right")

# ---------------------------------------------------------------------------
# Chart 4 — Absolute mapped counts
# ---------------------------------------------------------------------------

def plot_mapped_counts(ax: plt.Axes, df: pd.DataFrame) -> None:
    mean = df["Mapped_reads"].mean()
    bars = ax.bar(df["Label"], df["Mapped_reads"],
                  color=COLORS["mapped"], width=0.6, zorder=3)

    ax.axhline(mean, color=COLORS["threshold"],
               linewidth=1.5, linestyle="--", zorder=4)
    ax.text(len(df) - 0.5, mean * 1.01,
            f"Mean: {mean / 1e6:.1f}M",
            color=COLORS["threshold"], fontsize=9, ha="right")

    for bar, val in zip(bars, df["Mapped_reads"]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() * 1.01,
                f"{val / 1e6:.1f}M", ha="center", va="bottom", fontsize=8)

    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f"{v / 1e6:.0f}M")
    )
    ax.set_title("Absolute mapped read count per sample", fontsize=11,
                 fontweight="medium", pad=8)
    ax.set_ylabel("Mapped reads", fontsize=9)
    _style_ax(ax, df)

# ---------------------------------------------------------------------------
# Shared axis styling
# ---------------------------------------------------------------------------

def _style_ax(ax: plt.Axes, df: pd.DataFrame, set_xticks: bool = True) -> None:
    ax.set_xlim(-0.6, len(df) - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    if set_xticks:
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(df["Label"], rotation=30, ha="right", fontsize=9)
    ax.tick_params(axis="y", labelsize=9)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Plot sample QC metrics.")
    parser.add_argument("--input",     default="outputs/metrics.csv",  help="Path to metrics.csv")
    parser.add_argument("--output",    default="outputs/qc_plots.pdf",  help="Output file (.pdf, .png, .svg, …)")
    parser.add_argument("--threshold", default=0.60, type=float,
                        help="Low-coverage threshold as a fraction (default: 0.60)")
    parser.add_argument("--dpi",       default=150, type=int,   help="DPI for raster outputs")
    args = parser.parse_args()

    df = load_and_enrich(args.input, args.threshold)

    fig, axes = plt.subplots(4, 1, figsize=(10, 14))
    fig.suptitle("Sample QC summary", fontsize=13, fontweight="medium", y=0.98)
    plt.subplots_adjust(hspace=0.55)

    plot_mapping_rate(axes[0],       df, args.threshold * 100)
    plot_stacked_pct(axes[1],        df)
    plot_absolute_failures(axes[2],  df)
    plot_mapped_counts(axes[3], df)

    out = Path(args.output)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print(f"Saved → {out.resolve()}")
    plt.show()


if __name__ == "__main__":
    main()
