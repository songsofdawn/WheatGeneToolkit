import io
import math
import textwrap
from functools import lru_cache
from pathlib import Path
from typing import Set, Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def _get_file_cache_key(path: str) -> tuple[str, int, int]:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    return str(resolved), stat.st_mtime_ns, stat.st_size


@lru_cache(maxsize=16)
def _read_tsv_cached(path_str: str, _mtime_ns: int, _size: int) -> pd.DataFrame:
    return pd.read_csv(path_str, sep="\t", dtype=str)


def _read_tsv(path: str) -> pd.DataFrame:
    return _read_tsv_cached(*_get_file_cache_key(path)).copy()


def normalize_gene_list(gene_list: List[str]) -> List[str]:
    """去空、去重，保持原顺序"""
    cleaned = []
    seen = set()
    for gene in gene_list:
        g = str(gene).strip()
        if g and g not in seen:
            cleaned.append(g)
            seen.add(g)
    return cleaned


def load_term2gene(term2gene_path: str) -> pd.DataFrame:
    df = _read_tsv(term2gene_path)
    expected = {"go_id", "gene_id"}
    if not expected.issubset(df.columns):
        raise ValueError(f"{term2gene_path} 必须包含列: {expected}")
    return df[["go_id", "gene_id"]].dropna().drop_duplicates()


def load_term2name(term2name_path: str) -> pd.DataFrame:
    df = _read_tsv(term2name_path)
    expected = {"go_id", "go_term_name"}
    if not expected.issubset(df.columns):
        raise ValueError(f"{term2name_path} 必须包含列: {expected}")
    return df[["go_id", "go_term_name"]].dropna().drop_duplicates()


def load_go_metadata(metadata_path: str) -> pd.DataFrame:
    df = _read_tsv(metadata_path)
    expected = {"go_id", "ontology"}
    if not expected.issubset(df.columns):
        raise ValueError(f"{metadata_path} 必须包含列: {expected}")
    keep_cols = [c for c in ["go_id", "go_term_name", "go_namespace", "ontology"] if c in df.columns]
    return df[keep_cols].drop_duplicates()


def load_background_genes(background_path: str) -> Set[str]:
    df = _read_tsv(background_path)
    if "gene_id" not in df.columns:
        raise ValueError(f"{background_path} 必须包含列: gene_id")
    return set(df["gene_id"].dropna().astype(str))


def build_go_to_genes(term2gene: pd.DataFrame, background_genes: Set[str]) -> Dict[str, Set[str]]:
    go_to_genes: Dict[str, Set[str]] = {}
    for _, row in term2gene.iterrows():
        go_id = row["go_id"]
        gene_id = row["gene_id"]
        if gene_id in background_genes:
            go_to_genes.setdefault(go_id, set()).add(gene_id)
    return go_to_genes


def enrich_go(
    study_genes: Set[str],
    background_genes: Set[str],
    go_to_genes: Dict[str, Set[str]],
    min_geneset_size: int = 3,
    max_geneset_size: int = 2000,
) -> pd.DataFrame:
    study_genes = study_genes & background_genes

    N = len(background_genes)
    n = len(study_genes)

    if n == 0:
        raise ValueError("输入的 DEG 与背景基因集没有交集，无法进行富集分析。")

    results = []

    for go_id, term_genes in go_to_genes.items():
        M = len(term_genes)
        if M < min_geneset_size or M > max_geneset_size:
            continue

        overlap = study_genes & term_genes
        k = len(overlap)
        if k == 0:
            continue

        pvalue = hypergeom.sf(k - 1, N, M, n)

        results.append({
            "go_id": go_id,
            "BgRatio_numerator": M,
            "BgRatio_denominator": N,
            "GeneRatio_numerator": k,
            "GeneRatio_denominator": n,
            "Count": k,
            "pvalue": pvalue,
            "geneID": "/".join(sorted(overlap))
        })

    res_df = pd.DataFrame(results)

    if res_df.empty:
        return res_df

    reject, p_adjust, _, _ = multipletests(res_df["pvalue"], method="fdr_bh")
    res_df["p.adjust"] = p_adjust
    res_df["significant"] = reject

    res_df["GeneRatio"] = (
        res_df["GeneRatio_numerator"].astype(str) + "/" +
        res_df["GeneRatio_denominator"].astype(str)
    )
    res_df["BgRatio"] = (
        res_df["BgRatio_numerator"].astype(str) + "/" +
        res_df["BgRatio_denominator"].astype(str)
    )

    res_df["FoldEnrichment"] = (
        (res_df["GeneRatio_numerator"] / res_df["GeneRatio_denominator"]) /
        (res_df["BgRatio_numerator"] / res_df["BgRatio_denominator"])
    )

    res_df["minus_log10_padj"] = res_df["p.adjust"].apply(
        lambda x: -math.log10(x) if x > 0 else 300
    )
    res_df["minus_log10_pvalue"] = res_df["pvalue"].apply(
        lambda x: -math.log10(x) if x > 0 else 300
    )

    res_df = res_df.sort_values(
        ["p.adjust", "pvalue", "Count"],
        ascending=[True, True, False]
    ).reset_index(drop=True)

    return res_df


def finalize_result_table(
    res: pd.DataFrame,
    term2name: pd.DataFrame,
    metadata: pd.DataFrame
) -> pd.DataFrame:
    res = res.merge(term2name, on="go_id", how="left")
    res = res.merge(
        metadata.drop_duplicates(subset=["go_id"]),
        on="go_id",
        how="left",
        suffixes=("", "_meta")
    )

    if "go_term_name_meta" in res.columns:
        res["go_term_name"] = res["go_term_name"].fillna(res["go_term_name_meta"])
        res = res.drop(columns=["go_term_name_meta"])

    preferred_cols = [
        "go_id", "go_term_name", "ontology", "go_namespace",
        "Count", "GeneRatio", "BgRatio",
        "FoldEnrichment", "pvalue", "p.adjust",
        "minus_log10_pvalue", "minus_log10_padj",
        "geneID"
    ]
    other_cols = [c for c in res.columns if c not in preferred_cols]
    return res[preferred_cols + other_cols]


def build_summary_df(
    input_gene_count: int,
    background_gene_count: int,
    deg_in_background_count: int,
    tested_go_term_count: int,
    significant_go_term_count_fdr: int
) -> pd.DataFrame:
    return pd.DataFrame({
        "metric": [
            "input_deg_count",
            "background_gene_count",
            "deg_in_background_count",
            "tested_go_term_count",
            "significant_go_term_count_fdr"
        ],
        "value": [
            input_gene_count,
            background_gene_count,
            deg_in_background_count,
            tested_go_term_count,
            significant_go_term_count_fdr
        ]
    })


def wrap_go_terms(term, width: int = 35) -> str:
    """Wrap long GO labels without losing their meaning."""
    value = "" if pd.isna(term) else str(term)
    return "\n".join(textwrap.wrap(value, width=width, break_long_words=False)) or value


def _get_qvalue_col(df: pd.DataFrame, qvalue_col: str = "qvalue") -> Optional[str]:
    if qvalue_col in df.columns:
        return qvalue_col
    if "p.adjust" in df.columns:
        return "p.adjust"
    return None


def _scale_bubble_size(counts, min_size: int = 45, max_size: int = 360) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    if len(counts) == 0:
        return np.array([])
    if np.nanmax(counts) == np.nanmin(counts):
        return np.full(len(counts), (min_size + max_size) / 2)

    sqrt_counts = np.sqrt(counts)
    return min_size + (
        (sqrt_counts - sqrt_counts.min()) /
        (sqrt_counts.max() - sqrt_counts.min())
    ) * (max_size - min_size)


def _parse_ratio_to_float(value) -> float:
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if "/" not in text:
        return pd.to_numeric(text, errors="coerce")
    numerator, denominator = text.split("/", 1)
    numerator = pd.to_numeric(numerator, errors="coerce")
    denominator = pd.to_numeric(denominator, errors="coerce")
    if denominator == 0 or pd.isna(numerator) or pd.isna(denominator):
        return np.nan
    return float(numerator) / float(denominator)


def _prepare_go_plot_df(
    go_df: pd.DataFrame,
    ontology: str,
    top_n: int = 15,
    qvalue_col: str = "qvalue",
    count_col: str = "Count",
    label_wrap_width: int = 35,
) -> Tuple[pd.DataFrame, Optional[str], str]:
    if go_df is None or go_df.empty or "ontology" not in go_df.columns:
        return pd.DataFrame(), None, "Count"

    actual_qvalue_col = _get_qvalue_col(go_df, qvalue_col)
    if actual_qvalue_col is None or count_col not in go_df.columns:
        return pd.DataFrame(), actual_qvalue_col, "Count"

    plot_df = go_df[go_df["ontology"] == ontology].copy()
    if plot_df.empty:
        return plot_df, actual_qvalue_col, "Count"

    plot_df[actual_qvalue_col] = pd.to_numeric(plot_df[actual_qvalue_col], errors="coerce").fillna(1.0)
    plot_df[count_col] = pd.to_numeric(plot_df[count_col], errors="coerce").fillna(0)
    plot_df = plot_df[plot_df[count_col] > 0].copy()
    if plot_df.empty:
        return plot_df, actual_qvalue_col, "Count"

    plot_df = plot_df.sort_values(
        [actual_qvalue_col, count_col],
        ascending=[True, False],
    ).head(top_n)

    term_col = "go_term_name" if "go_term_name" in plot_df.columns else "go_id"
    plot_df["GO_term_wrapped"] = plot_df[term_col].fillna(plot_df.get("go_id", "")).apply(
        lambda value: wrap_go_terms(value, width=label_wrap_width)
    )

    x_col = "Count"
    if {"GeneRatio_numerator", "GeneRatio_denominator"}.issubset(plot_df.columns):
        denominator = pd.to_numeric(plot_df["GeneRatio_denominator"], errors="coerce").replace(0, np.nan)
        plot_df["GeneRatio_value"] = pd.to_numeric(
            plot_df["GeneRatio_numerator"], errors="coerce"
        ) / denominator
        if plot_df["GeneRatio_value"].notna().any():
            x_col = "GeneRatio_value"
    elif "GeneRatio" in plot_df.columns:
        plot_df["GeneRatio_value"] = plot_df["GeneRatio"].apply(_parse_ratio_to_float)
        if plot_df["GeneRatio_value"].notna().any():
            x_col = "GeneRatio_value"

    if "BgRatio_numerator" in plot_df.columns:
        bg_count = pd.to_numeric(plot_df["BgRatio_numerator"], errors="coerce").replace(0, np.nan)
        plot_df["RichFactor"] = plot_df[count_col] / bg_count
        if plot_df["RichFactor"].notna().any():
            x_col = "RichFactor"

    plot_df = plot_df.iloc[::-1].copy()
    return plot_df, actual_qvalue_col, x_col


def _qvalue_norm(values) -> plt.Normalize:
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if values.empty:
        return plt.Normalize(vmin=0, vmax=1)
    vmin = float(values.min())
    vmax = float(values.max())
    if vmax == vmin:
        vmin = 0.0
        vmax = max(vmax, 1e-6)
    return plt.Normalize(vmin=vmin, vmax=vmax)


def _add_count_legend(lax, counts, bubble_min_size: int, bubble_max_size: int):
    counts = pd.Series(counts).dropna().astype(int)
    if counts.empty:
        return

    legend_counts = sorted(set([
        int(counts.min()),
        int(counts.median()),
        int(counts.max()),
    ]))
    legend_sizes = _scale_bubble_size(
        legend_counts,
        min_size=bubble_min_size,
        max_size=bubble_max_size,
    )
    handles = [
        lax.scatter(
            [],
            [],
            s=size,
            facecolors="black",
            edgecolors="black",
            linewidths=0.45,
            alpha=0.95,
        )
        for size in legend_sizes
    ]
    lax.legend(
        handles,
        [str(count) for count in legend_counts],
        title="Count",
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.82, 0.5),
        scatterpoints=1,
        labelspacing=1.0,
        borderpad=0.2,
        handletextpad=0.8,
        fontsize=8,
        title_fontsize=9,
    )


def plot_go_barplot(
    go_df: pd.DataFrame,
    ontology: str,
    top_n: int = 15,
    qvalue_col: str = "qvalue",
    count_col: str = "Count",
    label_wrap_width: int = 35,
):
    plot_df, actual_qvalue_col, _ = _prepare_go_plot_df(
        go_df=go_df,
        ontology=ontology,
        top_n=top_n,
        qvalue_col=qvalue_col,
        count_col=count_col,
        label_wrap_width=label_wrap_width,
    )
    if plot_df.empty or actual_qvalue_col is None:
        return None

    height = max(6, 0.45 * len(plot_df))

    fig = plt.figure(figsize=(8.4, height))
    gs = fig.add_gridspec(
        nrows=5,
        ncols=2,
        width_ratios=[1.0, 0.055],
        height_ratios=[0.14, 0.46, 0.06, 0.20, 0.14],
        wspace=0.08,
        hspace=0.05,
    )

    ax = fig.add_subplot(gs[:, 0])
    cax = fig.add_subplot(gs[1, 1])

    norm = _qvalue_norm(plot_df[actual_qvalue_col])
    cmap = plt.cm.coolwarm_r

    ax.barh(
        plot_df["GO_term_wrapped"],
        plot_df[count_col],
        color=cmap(norm(plot_df[actual_qvalue_col].values)),
        edgecolor="none",
        height=0.75,
    )

    ax.set_title(f"{ontology} GO enrichment", fontsize=12, pad=10)
    ax.set_xlabel("Count", fontsize=12)
    ax.set_ylabel("")
    ax.set_facecolor("#EBEBEB")
    fig.patch.set_facecolor("white")
    ax.grid(True, axis="x", color="white", linewidth=1.2, alpha=0.95)
    ax.grid(True, axis="y", color="white", linewidth=0.8, alpha=0.85)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(float(plot_df[count_col].max()) * 1.12, 1))

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("qvalue", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    return fig


def plot_go_bubbleplot(
    go_df: pd.DataFrame,
    ontology: str,
    top_n: int = 15,
    qvalue_col: str = "qvalue",
    count_col: str = "Count",
    label_wrap_width: int = 35,
    bubble_min_size: int = 45,
    bubble_max_size: int = 360,
):
    plot_df, actual_qvalue_col, x_col = _prepare_go_plot_df(
        go_df=go_df,
        ontology=ontology,
        top_n=top_n,
        qvalue_col=qvalue_col,
        count_col=count_col,
        label_wrap_width=label_wrap_width,
    )
    if plot_df.empty or actual_qvalue_col is None:
        return None

    height = max(6, 0.45 * len(plot_df))

    fig = plt.figure(figsize=(8.4, height))
    gs = fig.add_gridspec(
        nrows=5,
        ncols=2,
        width_ratios=[1.0, 0.055],
        height_ratios=[0.14, 0.46, 0.06, 0.20, 0.14],
        wspace=0.08,
        hspace=0.05,
    )

    ax = fig.add_subplot(gs[:, 0])
    cax = fig.add_subplot(gs[1, 1])
    lax = fig.add_subplot(gs[3, 1])
    lax.axis("off")

    norm = _qvalue_norm(plot_df[actual_qvalue_col])
    cmap = plt.cm.coolwarm_r
    sizes = _scale_bubble_size(
        plot_df[count_col].values,
        min_size=bubble_min_size,
        max_size=bubble_max_size,
    )

    x_values = pd.to_numeric(plot_df[x_col], errors="coerce").fillna(plot_df[count_col])
    scatter = ax.scatter(
        x_values,
        plot_df["GO_term_wrapped"],
        s=sizes,
        c=plot_df[actual_qvalue_col],
        cmap=cmap,
        norm=norm,
        alpha=0.95,
        edgecolors="black",
        linewidths=0.45,
    )

    x_label = {
        "RichFactor": "RichFactor",
        "GeneRatio_value": "GeneRatio",
        "Count": "Count",
    }.get(x_col, x_col)

    ax.set_title(f"{ontology} GO bubble plot", fontsize=12, pad=10)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel("")
    ax.set_facecolor("#EBEBEB")
    fig.patch.set_facecolor("white")
    ax.grid(True, color="white", linewidth=1.2, alpha=0.95)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max(float(x_values.max()) * 1.14, 0.01))

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

    cbar = fig.colorbar(scatter, cax=cax)
    cbar.set_label("qvalue", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    _add_count_legend(lax, plot_df[count_col], bubble_min_size, bubble_max_size)

    return fig


def _figure_to_png_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def create_go_bubbleplot_bytes(
    df: pd.DataFrame,
    ontology: str,
    top_n: int = 15,
    qvalue_col: str = "qvalue",
    count_col: str = "Count",
    label_wrap_width: int = 35,
    bubble_min_size: int = 45,
    bubble_max_size: int = 360,
) -> Optional[bytes]:
    fig = plot_go_bubbleplot(
        go_df=df,
        ontology=ontology,
        top_n=top_n,
        qvalue_col=qvalue_col,
        count_col=count_col,
        label_wrap_width=label_wrap_width,
        bubble_min_size=bubble_min_size,
        bubble_max_size=bubble_max_size,
    )
    if fig is None:
        return None
    return _figure_to_png_bytes(fig)


def run_go_enrichment(
    gene_list: List[str],
    term2gene_path: str,
    term2name_path: str,
    metadata_path: str,
    background_path: str,
    min_size: int = 3,
    max_size: int = 2000,
    padj_cutoff: float = 0.05
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    gene_list = normalize_gene_list(gene_list)
    study_genes = set(gene_list)

    term2gene = load_term2gene(term2gene_path)
    term2name = load_term2name(term2name_path)
    metadata = load_go_metadata(metadata_path)
    background_genes = load_background_genes(background_path)

    overlap_input = study_genes & background_genes
    if len(overlap_input) == 0:
        raise ValueError("你的 DEG 列表与背景基因没有交集，请检查基因 ID 格式是否一致。")

    go_to_genes = build_go_to_genes(term2gene, background_genes)

    res = enrich_go(
        study_genes=study_genes,
        background_genes=background_genes,
        go_to_genes=go_to_genes,
        min_geneset_size=min_size,
        max_geneset_size=max_size
    )

    if res.empty:
        summary_df = build_summary_df(
            input_gene_count=len(study_genes),
            background_gene_count=len(background_genes),
            deg_in_background_count=len(overlap_input),
            tested_go_term_count=0,
            significant_go_term_count_fdr=0
        )
        return res, res.copy(), summary_df

    res = finalize_result_table(res, term2name, metadata)
    sig = res[res["p.adjust"] <= padj_cutoff].copy()

    summary_df = build_summary_df(
        input_gene_count=len(study_genes),
        background_gene_count=len(background_genes),
        deg_in_background_count=len(overlap_input),
        tested_go_term_count=len(res),
        significant_go_term_count_fdr=(res["p.adjust"] <= padj_cutoff).sum()
    )

    return res, sig, summary_df


def create_go_barplot_bytes(
    df: pd.DataFrame,
    top_n: int = 15,
    plot_metric: str = "qvalue",
    figsize=None
) -> Optional[bytes]:
    if df.empty:
        return None

    df = df[df["ontology"].isin(["BP", "CC", "MF"])].copy()
    if df.empty:
        return None

    # GO 图只使用 qvalue 表示显著性；当前结果表中 qvalue 等价于 p.adjust。
    qvalue_col = "qvalue" if "qvalue" in df.columns else "p.adjust"
    if qvalue_col not in df.columns:
        return None

    subsets = []
    for onto in ["BP", "CC", "MF"]:
        sub = df[df["ontology"] == onto].copy()
        if sub.empty:
            continue
        sub[qvalue_col] = pd.to_numeric(sub[qvalue_col], errors="coerce").fillna(1.0)
        sub["Count"] = pd.to_numeric(sub["Count"], errors="coerce").fillna(0)
        sub = sub.sort_values(qvalue_col, ascending=True).head(top_n)
        sub["Description_short"] = sub["go_term_name"].fillna(sub["go_id"]).apply(
            lambda value: "\n".join(textwrap.wrap(str(value), width=42))
        )
        subsets.append((onto, sub))

    if not subsets:
        return None

    n_panels = len(subsets)
    total_terms = sum(len(sub) for _, sub in subsets)
    fig_height = max(5.5, 0.42 * total_terms + 1.2 * n_panels)
    fig, axes = plt.subplots(
        n_panels,
        1,
        figsize=figsize or (9.5, fig_height),
        constrained_layout=False,
    )

    if n_panels == 1:
        axes = [axes]

    all_qvalues = pd.concat([sub[qvalue_col] for _, sub in subsets], ignore_index=True)
    vmin = float(all_qvalues.min()) if not all_qvalues.empty else 0.0
    vmax = float(all_qvalues.max()) if not all_qvalues.empty else 1.0
    if vmax == vmin:
        vmin = 0.0
        vmax = max(vmax, 1e-6)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.cm.coolwarm_r

    for ax, (onto, sub) in zip(axes, subsets):
        sub = sub.sort_values(qvalue_col, ascending=False)
        ax.barh(
            sub["Description_short"],
            sub["Count"],
            color=cmap(norm(sub[qvalue_col].values)),
            edgecolor="none",
        )
        ax.set_title(f"{onto} GO enrichment", fontsize=12)
        ax.set_xlabel("Count", fontsize=11)
        ax.set_ylabel("GO term")
        ax.set_facecolor("#EBEBEB")
        ax.grid(True, axis="x", color="white", linewidth=1.2)
        ax.grid(True, axis="y", color="white", linewidth=0.8, alpha=0.85)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(axis="both", labelsize=9)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.subplots_adjust(left=0.34, right=0.82, top=0.95, bottom=0.08, hspace=0.48)
    cax = fig.add_axes([0.86, 0.42, 0.024, 0.22])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("qvalue", fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()


def create_go_barplot_panel_bytes(
    df: pd.DataFrame,
    ontology: str,
    top_n: int = 15,
    qvalue_col: str = "qvalue",
    count_col: str = "Count",
    label_wrap_width: int = 35,
) -> Optional[bytes]:
    fig = plot_go_barplot(
        go_df=df,
        ontology=ontology,
        top_n=top_n,
        qvalue_col=qvalue_col,
        count_col=count_col,
        label_wrap_width=label_wrap_width,
    )
    if fig is None:
        return None
    return _figure_to_png_bytes(fig)
