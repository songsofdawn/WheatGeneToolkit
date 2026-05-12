import io
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


NA_VALUES = ["NA", "NaN", "nan", ""]


def _read_text_from_input(file_or_text) -> str:
    if file_or_text is None:
        return ""

    if hasattr(file_or_text, "read"):
        data = file_or_text.read()
        if isinstance(data, bytes):
            return data.decode("utf-8-sig", errors="replace")
        return str(data)

    if isinstance(file_or_text, bytes):
        return file_or_text.decode("utf-8-sig", errors="replace")

    return str(file_or_text)


def _detect_sep(text: str, sep: Optional[str] = None):
    if sep is not None:
        return sep

    first_data_line = ""
    for line in text.splitlines():
        if line.strip():
            first_data_line = line
            break

    if "\t" in first_data_line:
        return "\t"
    if "," in first_data_line:
        return ","
    return r"\s+"


def read_volcano_table(file_or_text, sep=None) -> pd.DataFrame:
    text = _read_text_from_input(file_or_text)
    if not text.strip():
        return pd.DataFrame()

    detected_sep = _detect_sep(text, sep=sep)
    return pd.read_csv(
        io.StringIO(text),
        sep=detected_sep,
        engine="python",
        header=0,
        na_values=NA_VALUES,
        keep_default_na=True,
    )


def _first_exact_match(columns, candidates):
    lower_map = {str(col).lower(): col for col in columns}
    for candidate in candidates:
        col = lower_map.get(candidate.lower())
        if col is not None:
            return col
    return None


def _first_contains_match(columns, patterns):
    lowered = [(col, str(col).lower()) for col in columns]
    for pattern in patterns:
        p = pattern.lower()
        for col, lower_col in lowered:
            if p in lower_col:
                return col
    return None


def auto_detect_columns(df) -> dict:
    columns = list(df.columns)

    gene_id_col = _first_exact_match(
        columns,
        ["st_gene_id", "gene_id", "GeneID", "gene", "Gene", "id", "ID"],
    )

    log2fc_col = _first_contains_match(
        columns,
        [
            "diffexp_log2fc",
            "log2foldchange",
            "log2_fc",
            "log2fc",
            "fold_change",
            "foldchange",
        ],
    )

    pvalue_col = _first_contains_match(
        columns,
        [
            "deseq2_pvalue",
            "pvalue",
            "p_value",
            "p.value",
            "pval",
            "padj",
            "qvalue",
            "fdr",
        ],
    )

    return {
        "gene_id_col": gene_id_col,
        "log2fc_col": log2fc_col,
        "pvalue_col": pvalue_col,
    }


def infer_auto_y_cap(values, hard_max=320.0):
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if values.empty:
        return None

    max_value = float(values.max())
    if max_value <= hard_max:
        return None
    return hard_max


def _auto_soft_y_limit(values, quantile=0.995, min_limit=10.0):
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if values.empty:
        return min_limit
    return max(float(values.quantile(quantile)), min_limit)


def _soft_compress_positive(values, linear_limit):
    values = np.asarray(values, dtype=float)
    limit = max(float(linear_limit), 1e-9)
    compressed = values.copy()
    mask = values > limit
    compressed[mask] = limit + np.log1p(values[mask] - limit)
    return compressed


def _soft_compress_symmetric(values, linear_limit):
    values = np.asarray(values, dtype=float)
    limit = max(float(linear_limit), 1e-9)
    signs = np.sign(values)
    abs_values = np.abs(values)
    compressed = abs_values.copy()
    mask = abs_values > limit
    compressed[mask] = limit + np.log1p(abs_values[mask] - limit)
    return signs * compressed


def _auto_soft_x_limit(values, quantile=0.995, min_limit=2.0):
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if values.empty:
        return min_limit
    return max(float(values.abs().quantile(quantile)), min_limit)


def prepare_volcano_data(
    df,
    gene_id_col=None,
    log2fc_col=None,
    pvalue_col=None,
    log2fc_cutoff=1.0,
    pvalue_cutoff=0.05,
    use_abs_log2fc=True,
    max_neg_log10_pvalue=None,
):
    if df is None or df.empty:
        return pd.DataFrame()

    if log2fc_col is None or pvalue_col is None:
        raise ValueError("log2fc_col and pvalue_col are required")

    if gene_id_col is not None and gene_id_col in df.columns:
        gene_ids = df[gene_id_col].astype(str).str.strip()
    else:
        gene_ids = pd.Series([f"row_{i + 1}" for i in range(len(df))], index=df.index)

    result = pd.DataFrame({
        "gene_id": gene_ids,
        "log2FC": pd.to_numeric(df[log2fc_col], errors="coerce"),
        "pvalue": pd.to_numeric(df[pvalue_col], errors="coerce"),
    })

    result = result.replace({"gene_id": {"": np.nan, "nan": np.nan, "NA": np.nan}})
    result = result.dropna(subset=["gene_id", "log2FC", "pvalue"]).copy()
    result = result[(result["pvalue"] > 0) & (result["pvalue"] <= 1)].copy()

    if result.empty:
        return result

    result["neg_log10_pvalue_raw"] = -np.log10(result["pvalue"].clip(lower=1e-300))

    if max_neg_log10_pvalue is None:
        max_neg_log10_pvalue = infer_auto_y_cap(result["neg_log10_pvalue_raw"])

    if max_neg_log10_pvalue is not None:
        cap = float(max_neg_log10_pvalue)
        result["neg_log10_pvalue_plot"] = result["neg_log10_pvalue_raw"].clip(upper=cap)
        result["is_y_capped"] = result["neg_log10_pvalue_raw"] > cap
    else:
        result["neg_log10_pvalue_plot"] = result["neg_log10_pvalue_raw"]
        result["is_y_capped"] = False

    soft_y_limit = _auto_soft_y_limit(result["neg_log10_pvalue_plot"])
    result["neg_log10_pvalue_display"] = _soft_compress_positive(
        result["neg_log10_pvalue_plot"],
        soft_y_limit,
    )
    result["is_y_compressed"] = result["neg_log10_pvalue_plot"] > soft_y_limit

    if use_abs_log2fc:
        up_mask = (result["log2FC"] >= log2fc_cutoff) & (result["pvalue"] <= pvalue_cutoff)
        down_mask = (result["log2FC"] <= -log2fc_cutoff) & (result["pvalue"] <= pvalue_cutoff)
    else:
        up_mask = (result["log2FC"] >= log2fc_cutoff) & (result["pvalue"] <= pvalue_cutoff)
        down_mask = (result["log2FC"] <= -log2fc_cutoff) & (result["pvalue"] <= pvalue_cutoff)

    result["regulation"] = "Not significant"
    result.loc[up_mask, "regulation"] = "Up-regulated"
    result.loc[down_mask, "regulation"] = "Down-regulated"

    return result.reset_index(drop=True)


def _get_xlim(values):
    values = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if values.empty:
        return (-2, 2)

    max_abs = max(float(values.abs().max()), 1.0)
    abs_limit = max_abs * 1.08
    return (-abs_limit, abs_limit)


def _make_symmetric_ticks(raw_limit, display_limit, soft_limit):
    candidates = [-raw_limit, -soft_limit, -1.0, 0.0, 1.0, soft_limit, raw_limit]
    ticks = []
    labels = []
    for value in candidates:
        if abs(value) > raw_limit + 1e-9:
            continue
        if any(abs(value - existing) < 1e-9 for existing in ticks):
            continue
        ticks.append(value)
        labels.append(f"{value:g}")
    tick_positions = _soft_compress_symmetric(np.array(ticks), soft_limit)
    keep = np.abs(tick_positions) <= display_limit + 1e-9
    return tick_positions[keep], [label for label, ok in zip(labels, keep) if ok]


def _make_positive_ticks(raw_max, display_max, soft_limit):
    candidates = [0.0, 1.0, 2.0, 5.0, soft_limit, raw_max]
    ticks = []
    labels = []
    for value in candidates:
        if value > raw_max + 1e-9:
            continue
        if any(abs(value - existing) < 1e-9 for existing in ticks):
            continue
        ticks.append(value)
        labels.append(f"{value:g}")
    tick_positions = _soft_compress_positive(np.array(ticks), soft_limit)
    keep = tick_positions <= display_max + 1e-9
    return tick_positions[keep], [label for label, ok in zip(labels, keep) if ok]


def _label_top_genes(ax, volcano_df, top_label_n):
    if top_label_n <= 0 or volcano_df.empty:
        return

    label_df = volcano_df.copy()
    label_df["rank_score"] = label_df["neg_log10_pvalue_raw"] * label_df["log2FC"].abs()
    label_df = label_df.sort_values(["rank_score", "neg_log10_pvalue_raw"], ascending=False).head(top_label_n)

    for idx, (_, row) in enumerate(label_df.iterrows()):
        x = float(row["log2FC_display"] if "log2FC_display" in row else row["log2FC"])
        y = float(row["neg_log10_pvalue_display"] if "neg_log10_pvalue_display" in row else row["neg_log10_pvalue_plot"])
        align = "left" if x >= 0 else "right"
        dx = 5 if x >= 0 else -5
        dy = 4 + (idx % 3) * 3
        ax.annotate(
            str(row["gene_id"]),
            xy=(x, y),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=align,
            va="bottom",
            fontsize=7.5,
            color="#333333",
            arrowprops={"arrowstyle": "-", "color": "#9A9A9A", "linewidth": 0.5},
        )


def plot_volcano(
    volcano_df,
    log2fc_cutoff=1.0,
    pvalue_cutoff=0.05,
    title="Volcano Plot",
    top_label_n=10,
    label_mode="top_significant",
    max_neg_log10_pvalue=None,
    point_size=14,
    alpha=0.75,
):
    if volcano_df is None or volcano_df.empty:
        raise ValueError("volcano_df is empty")

    plot_df = volcano_df.copy()
    soft_x_limit = _auto_soft_x_limit(plot_df["log2FC"])
    soft_y_limit = _auto_soft_y_limit(plot_df["neg_log10_pvalue_plot"])
    plot_df["log2FC_display"] = _soft_compress_symmetric(plot_df["log2FC"], soft_x_limit)
    plot_df["neg_log10_pvalue_display"] = _soft_compress_positive(plot_df["neg_log10_pvalue_plot"], soft_y_limit)
    plot_df["is_x_compressed"] = plot_df["log2FC"].abs() > soft_x_limit
    plot_df["is_y_compressed"] = plot_df["neg_log10_pvalue_plot"] > soft_y_limit
    log2fc_cutoff_display = float(_soft_compress_symmetric(np.array([log2fc_cutoff]), soft_x_limit)[0])
    neg_log2fc_cutoff_display = float(_soft_compress_symmetric(np.array([-log2fc_cutoff]), soft_x_limit)[0])
    y_cutoff_raw = -np.log10(max(float(pvalue_cutoff), 1e-300))
    y_cutoff_display = float(_soft_compress_positive(np.array([y_cutoff_raw]), soft_y_limit)[0])

    fig, ax = plt.subplots(figsize=(8.8, 7.0))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#FAFAFA")

    styles = {
        "Not significant": {"color": "#BDBDBD", "label": "NotSig", "zorder": 1, "alpha": min(alpha, 0.42)},
        "Down-regulated": {"color": "#4C78A8", "label": "Down", "zorder": 2, "alpha": alpha},
        "Up-regulated": {"color": "#E53935", "label": "Up", "zorder": 3, "alpha": alpha},
    }

    for regulation in ["Not significant", "Down-regulated", "Up-regulated"]:
        sub = plot_df[plot_df["regulation"] == regulation].copy()
        if sub.empty:
            continue
        style = styles[regulation]
        ax.scatter(
            sub["log2FC_display"],
            sub["neg_log10_pvalue_display"],
            s=point_size,
            c=style["color"],
            alpha=style["alpha"],
            edgecolors="none",
            label=style["label"],
            zorder=style["zorder"],
            rasterized=len(plot_df) > 30000,
        )

    capped = plot_df[plot_df["is_y_capped"]].copy()
    if not capped.empty:
        cap_y = float(capped["neg_log10_pvalue_display"].max())
        ax.axhline(cap_y, color="#A0A0A0", linestyle=":", linewidth=0.8, alpha=0.75)

    ax.axvline(log2fc_cutoff_display, color="#222222", linestyle=(0, (4, 4)), linewidth=1.15, alpha=0.9)
    ax.axvline(neg_log2fc_cutoff_display, color="#222222", linestyle=(0, (4, 4)), linewidth=1.15, alpha=0.9)
    ax.axhline(y_cutoff_display, color="#222222", linestyle=(0, (4, 4)), linewidth=1.15, alpha=0.9)

    up_count = int((plot_df["regulation"] == "Up-regulated").sum())
    down_count = int((plot_df["regulation"] == "Down-regulated").sum())
    sig_count = up_count + down_count
    subtitle = f"Significant DEG = {sig_count} | Up = {up_count} | Down = {down_count}"
    ax.set_title(f"{title}\n{subtitle}", fontsize=15, loc="left", pad=12)
    ax.set_xlabel(r"$\log_2(\mathrm{FoldChange})$", fontsize=12)
    ax.set_ylabel(r"$-\log_{10}(p\mathrm{-value})$", fontsize=12)
    ax.grid(True, which="major", color="#E6E6E6", linewidth=0.9, alpha=0.95)
    ax.minorticks_on()
    ax.grid(True, which="minor", color="#F0F0F0", linewidth=0.55, alpha=0.9)
    ax.set_axisbelow(True)

    for spine in ax.spines.values():
        spine.set_color("#4D4D4D")
        spine.set_linewidth(0.9)

    x_raw_limit = max(float(plot_df["log2FC"].abs().max()), abs(log2fc_cutoff) * 1.4, 1.0)
    x_display_limit = float(_soft_compress_symmetric(np.array([x_raw_limit]), soft_x_limit)[0]) * 1.06
    ax.set_xlim(-x_display_limit, x_display_limit)
    x_ticks, x_labels = _make_symmetric_ticks(x_raw_limit, x_display_limit, soft_x_limit)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)

    y_raw_max = max(float(plot_df["neg_log10_pvalue_plot"].max()), y_cutoff_raw * 1.25, 1.0)
    y_display_max = float(_soft_compress_positive(np.array([y_raw_max]), soft_y_limit)[0]) * 1.08
    ax.set_ylim(0, y_display_max)
    y_ticks, y_labels = _make_positive_ticks(y_raw_max, y_display_max, soft_y_limit)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)

    if label_mode == "top_significant":
        _label_top_genes(ax, plot_df, top_label_n)

    if not capped.empty:
        ax.text(
            0.99,
            0.98,
            f"{len(capped)} capped at top",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            color="#555555",
        )

    x_compressed_count = int(plot_df["is_x_compressed"].sum())
    y_compressed_count = int(plot_df["is_y_compressed"].sum())
    if x_compressed_count or y_compressed_count:
        ax.text(
            0.01,
            0.98,
            f"axis compressed: x={x_compressed_count}, y={y_compressed_count}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="#555555",
        )

    ax.legend(
        title="Group",
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=10,
        title_fontsize=11,
        markerscale=1.2,
    )

    fig.subplots_adjust(left=0.12, right=0.80, top=0.88, bottom=0.12)
    return fig


def figure_to_png_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()
