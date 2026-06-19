from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ROOT = Path(r"D:\TP53_project")
RESULTS_FIG = ROOT / "results" / "figures"


def _require(pkg: str) -> None:
    __import__(pkg)


def _savefig(fig, out_path: Path, dpi: int = 300) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")


def figure1_lollipop_hotspots() -> Path:
    _require("matplotlib")
    import matplotlib.pyplot as plt

    df = pd.read_csv(ROOT / "data" / "mutations" / "tp53_hotspot_frequency.csv")
    df = df.sort_values(["Mutation_Frequency", "Position"], ascending=[False, True]).reset_index(drop=True)
    top = df.head(30).copy()

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.vlines(top["Position"], [0], top["Mutation_Frequency"], color="#2c3e50", linewidth=1.8, alpha=0.9)
    ax.scatter(top["Position"], top["Mutation_Frequency"], s=65, color="#c0392b", zorder=3)

    ax.set_title("TP53 hotspot positions (top 30; ClinVar filtered SNVs)", fontsize=13)
    ax.set_xlabel("Amino acid position (TP53)")
    ax.set_ylabel("Frequency in filtered ClinVar dataset")
    ax.grid(axis="y", linestyle="--", alpha=0.25)

    # Annotate the top hotspot (expected to be 281 in this project)
    if not top.empty:
        best = top.iloc[0]
        ax.annotate(
            f"Top hotspot: {int(best['Position'])} (n={int(best['Mutation_Frequency'])})",
            xy=(best["Position"], best["Mutation_Frequency"]),
            xytext=(best["Position"] + 10, best["Mutation_Frequency"] + 1),
            arrowprops=dict(arrowstyle="->", lw=1),
            fontsize=10,
        )

    out = RESULTS_FIG / "Figure1_TP53_hotspot_lollipop.png"
    _savefig(fig, out)
    plt.close(fig)
    return out


def figure2_volcano() -> Path:
    _require("matplotlib")
    import matplotlib.pyplot as plt

    res_path = ROOT / "data" / "expression" / "DEA_results" / "DESeq2_all_results_annotated.csv"
    df = pd.read_csv(res_path)

    # Expected columns from DESeq2 CSV: log2FoldChange, padj, SYMBOL/ENSEMBL
    df = df.dropna(subset=["log2FoldChange", "padj"]).copy()
    df["neglog10padj"] = -df["padj"].apply(lambda p: math.log10(p) if p > 0 else float("inf"))

    sig = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() > 1)

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    ax.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "neglog10padj"], s=6, c="#95a5a6", alpha=0.35, linewidths=0)
    ax.scatter(df.loc[sig, "log2FoldChange"], df.loc[sig, "neglog10padj"], s=7, c="#e74c3c", alpha=0.55, linewidths=0)

    ax.axvline(1, color="#34495e", lw=1, linestyle="--", alpha=0.7)
    ax.axvline(-1, color="#34495e", lw=1, linestyle="--", alpha=0.7)
    ax.axhline(-math.log10(0.05), color="#34495e", lw=1, linestyle="--", alpha=0.7)

    ax.set_title("Volcano plot: TP53-mutant vs TP53-wildtype (TCGA-BRCA)", fontsize=12.5)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.grid(alpha=0.15)

    out = RESULTS_FIG / "Figure2_Volcano_TP53_mutant_vs_wildtype.png"
    _savefig(fig, out)
    plt.close(fig)
    return out


def figure3_heatmap_group_means() -> Path:
    _require("matplotlib")
    _require("seaborn")
    import matplotlib.pyplot as plt
    import seaborn as sns

    sig_path = ROOT / "data" / "expression" / "DEA_results" / "DESeq2_significant_DEGs_annotated.csv"
    norm_path = ROOT / "data" / "expression" / "DEA_results" / "normalized_counts.csv"
    meta_path = ROOT / "data" / "expression" / "deseq2_ready_metadata.csv"

    sig = pd.read_csv(sig_path)
    sig = sig.dropna(subset=["ENSEMBL", "padj"]).sort_values("padj", ascending=True)
    top50 = sig.head(50)["ENSEMBL"].astype(str).tolist()

    norm = pd.read_csv(norm_path, index_col=0)
    meta = pd.read_csv(meta_path)

    # metadata columns: sample_barcode, case_id, TP53_status
    # normalized_counts columns are sample barcodes
    mutant_samples = meta.loc[meta["TP53_status"] == "TP53_mutant", "sample_barcode"].astype(str).tolist()
    wt_samples = meta.loc[meta["TP53_status"] == "TP53_wildtype", "sample_barcode"].astype(str).tolist()

    mutant_samples = [s for s in mutant_samples if s in norm.columns]
    wt_samples = [s for s in wt_samples if s in norm.columns]

    mat = norm.loc[norm.index.intersection(top50), mutant_samples + wt_samples].copy()

    # Group means to keep the plot readable and deterministic across N~1000 samples
    mean_mut = mat[mutant_samples].mean(axis=1)
    mean_wt = mat[wt_samples].mean(axis=1)
    group_mat = pd.DataFrame({"TP53_mutant_mean": mean_mut, "TP53_wildtype_mean": mean_wt})

    # z-score per gene
    group_mat = group_mat.sub(group_mat.mean(axis=1), axis=0)
    group_mat = group_mat.div(group_mat.std(axis=1).replace(0, pd.NA), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(6.4, 10.5))
    sns.heatmap(group_mat, cmap="vlag", center=0, ax=ax, cbar_kws={"label": "Z-score (per gene)"})
    ax.set_title("Top 50 DEGs heatmap (group means; z-scored per gene)", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("Gene (Ensembl ID)")

    out = RESULTS_FIG / "Figure3_Heatmap_Top50_DEGs_group_means.png"
    _savefig(fig, out)
    plt.close(fig)
    return out


def figure4_hub_top20_svg() -> Path:
    _require("matplotlib")
    import matplotlib.pyplot as plt

    hub_path = ROOT / "data" / "network_analysis" / "ranked_HUB.csv"
    # File contains a descriptive first line before the CSV header.
    df = pd.read_csv(hub_path, skiprows=1)
    df = df.sort_values(["Score", "Name"], ascending=[False, True]).head(20).copy()
    df = df.iloc[::-1]  # for horizontal bar plot (top at bottom -> reversed)

    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    ax.barh(df["Name"], df["Score"], color="#2c3e50", alpha=0.9)
    ax.set_title("Top 20 hub genes ranked by MCC (CytoHubba)", fontsize=12.5)
    ax.set_xlabel("MCC score")
    ax.set_ylabel("")
    ax.grid(axis="x", linestyle="--", alpha=0.25)

    out = RESULTS_FIG / "Figure4_HUB_TOP20_MCC.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, format="svg", bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    RESULTS_FIG.mkdir(parents=True, exist_ok=True)

    outs = [
        figure1_lollipop_hotspots(),
        figure2_volcano(),
        figure3_heatmap_group_means(),
        figure4_hub_top20_svg(),
    ]

    # Copy existing structural SVGs into results/figures for consistent thesis paths
    struct_dir = ROOT / "data" / "structure_prediction"
    for name in [
        "psipredChart_WT.svg",
        "psipredChart_Mutant.svg",
        "annotationGrid_WT.svg",
        "annotationGrid_Mutant.svg",
    ]:
        src = struct_dir / name
        if src.exists():
            dst = RESULTS_FIG / name
            dst.write_bytes(src.read_bytes())
            outs.append(dst)

    print("Generated/collected thesis figures:")
    for p in outs:
        print(f"- {p}")


if __name__ == "__main__":
    main()

