# -*- coding: utf-8 -*-

"""
KEGG 富集分析工具函数

适配目录：
data/kegg_mapping/
├── gene2ko_clean.tsv
├── kegg_ko2pathway.tsv
└── kegg_pathway2name.tsv

核心功能：
1. 读取 gene -> KO 映射
2. 读取 KO -> KEGG pathway 映射
3. 使用 KO 层面的超几何检验进行 KEGG 富集
4. 输出全部结果、显著结果、分析摘要
5. 生成 KEGG 气泡图和条形图的 PNG bytes，供 Streamlit 显示和下载
"""

import io
import textwrap
from functools import lru_cache
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import hypergeom


def _get_file_cache_key(path: str) -> tuple[str, int, int]:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    return str(resolved), stat.st_mtime_ns, stat.st_size


@lru_cache(maxsize=16)
def _read_tsv_cached(path_str: str, _mtime_ns: int, _size: int) -> pd.DataFrame:
    return pd.read_csv(path_str, sep="\t", dtype=str)


def _read_tsv(path: str) -> pd.DataFrame:
    return _read_tsv_cached(*_get_file_cache_key(path)).copy()


# ==============================
# 基础工具函数
# ==============================

def bh_adjust(pvalues):
    """
    Benjamini-Hochberg FDR 校正
    """
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)

    if n == 0:
        return np.array([])

    order = np.argsort(pvalues)
    ranked = pvalues[order]

    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    result = np.empty(n)
    result[order] = adjusted

    return result


def wrap_text(text, width=42):
    """
    通路名称自动换行
    """
    return "\n".join(textwrap.wrap(str(text), width=width))


def scale_bubble_size(counts, min_size=25, max_size=220):
    """
    对 Count 做 sqrt 缩放，避免大通路气泡过大
    """
    counts = np.asarray(counts, dtype=float)

    if len(counts) == 0:
        return np.array([])

    if counts.max() == counts.min():
        return np.full(len(counts), (min_size + max_size) / 2)

    sqrt_counts = np.sqrt(counts)

    sizes = min_size + (
        (sqrt_counts - sqrt_counts.min()) /
        (sqrt_counts.max() - sqrt_counts.min())
    ) * (max_size - min_size)

    return sizes


def get_x_cap(values, clip=True, mode="quantile", quantile=0.95, fixed_value=30, min_x_cap=5):
    """
    获取 -log₁₀(P-value) 的截断上限。
    截断是为了避免极端显著通路把横轴拉得过长。
    """
    values = np.asarray(values, dtype=float)

    if len(values) == 0:
        return 1

    if not clip:
        return values.max()

    if mode == "fixed":
        cap = fixed_value
    else:
        cap = np.quantile(values, quantile)

    cap = max(cap, min_x_cap)

    return cap


# ==============================
# 数据读取
# ==============================

def load_gene2ko(gene2ko_path):
    """
    读取 gene2ko_clean.tsv，并处理一个基因对应多个 KO 的情况。

    输入文件要求：
    gene_id    ko_list
    Traes...   K00001
    Traes...   K13343;K16284

    返回：
    background DataFrame:
    gene       KO
    Traes...   K00001
    Traes...   K13343
    Traes...   K16284
    """

    df = _read_tsv(gene2ko_path)

    required_cols = {"gene_id", "ko_list"}
    if not required_cols.issubset(set(df.columns)):
        raise ValueError(
            f"{gene2ko_path} 必须包含列 gene_id 和 ko_list，当前列名为：{list(df.columns)}"
        )

    df = df[["gene_id", "ko_list"]].copy()
    df.columns = ["gene", "KO"]

    df = df.dropna()
    df = df[(df["gene"] != "") & (df["KO"] != "")]

    # 拆分多个 KO
    df["KO"] = df["KO"].str.split(";")
    df = df.explode("KO")

    df["gene"] = df["gene"].astype(str).str.strip()
    df["KO"] = df["KO"].astype(str).str.strip()

    df = df[(df["gene"] != "") & (df["KO"] != "")]
    df = df[df["KO"].str.startswith("K")]

    df = df.drop_duplicates()

    return df


def load_kegg_mapping(ko2pathway_path, pathway2name_path):
    """
    读取本地 KEGG 注释文件。

    kegg_ko2pathway.tsv:
        KO      pathway
        K00001  ko00010

    kegg_pathway2name.tsv:
        pathway     Description
        ko00010     Glycolysis / Gluconeogenesis
    """

    ko2pathway = _read_tsv(ko2pathway_path)
    pathway2name = _read_tsv(pathway2name_path)

    required_cols_1 = {"KO", "pathway"}
    required_cols_2 = {"pathway", "Description"}

    if not required_cols_1.issubset(set(ko2pathway.columns)):
        raise ValueError(
            f"{ko2pathway_path} 必须包含列 KO 和 pathway，当前列名为：{list(ko2pathway.columns)}"
        )

    if not required_cols_2.issubset(set(pathway2name.columns)):
        raise ValueError(
            f"{pathway2name_path} 必须包含列 pathway 和 Description，当前列名为：{list(pathway2name.columns)}"
        )

    ko2pathway = ko2pathway.dropna()
    pathway2name = pathway2name.dropna()

    ko2pathway["KO"] = ko2pathway["KO"].astype(str).str.strip()
    ko2pathway["pathway"] = ko2pathway["pathway"].astype(str).str.strip()

    pathway2name["pathway"] = pathway2name["pathway"].astype(str).str.strip()
    pathway2name["Description"] = pathway2name["Description"].astype(str).str.strip()

    ko2pathway = ko2pathway[
        ko2pathway["KO"].str.startswith("K") &
        ko2pathway["pathway"].str.startswith("ko")
    ]

    ko2pathway = ko2pathway.drop_duplicates()

    pathway_to_kos = (
        ko2pathway
        .groupby("pathway")["KO"]
        .apply(lambda x: set(x))
        .to_dict()
    )

    pathway_to_name = dict(
        zip(pathway2name["pathway"], pathway2name["Description"])
    )

    return pathway_to_kos, pathway_to_name, ko2pathway, pathway2name


# ==============================
# KEGG 富集分析
# ==============================

def add_gene_name_columns(result_df, background, deg_genes):
    """
    给富集结果添加：
    AllGeneNames
    DEG_GeneNames
    """

    if result_df.empty:
        return result_df

    ko_to_all_genes = (
        background
        .groupby("KO")["gene"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    deg_background = background[background["gene"].isin(deg_genes)].copy()

    ko_to_deg_genes = (
        deg_background
        .groupby("KO")["gene"]
        .apply(lambda x: sorted(set(x)))
        .to_dict()
    )

    all_gene_names_list = []
    deg_gene_names_list = []

    for _, row in result_df.iterrows():
        kos_in_pathway = str(row["geneID"]).split("/")

        all_genes = []
        deg_genes_in_pathway = []

        for ko in kos_in_pathway:
            all_genes.extend(ko_to_all_genes.get(ko, []))
            deg_genes_in_pathway.extend(ko_to_deg_genes.get(ko, []))

        all_gene_names_list.append("; ".join(sorted(set(all_genes))))
        deg_gene_names_list.append("; ".join(sorted(set(deg_genes_in_pathway))))

    result_df = result_df.copy()
    result_df["AllGeneNames"] = all_gene_names_list
    result_df["DEG_GeneNames"] = deg_gene_names_list

    return result_df


def run_kegg_enrichment(
    gene_list,
    gene2ko_path,
    ko2pathway_path,
    pathway2name_path,
    min_size=3,
    max_size=500,
    pvalue_cutoff=0.05
):
    """
    KEGG 富集主函数。

    参数：
    gene_list:
        DEG 基因列表，一行一个基因号

    gene2ko_path:
        data/kegg_mapping/gene2ko_clean.tsv

    ko2pathway_path:
        data/kegg_mapping/kegg_ko2pathway.tsv

    pathway2name_path:
        data/kegg_mapping/kegg_pathway2name.tsv

    min_size:
        pathway 最小 KO 数

    max_size:
        pathway 最大 KO 数

    pvalue_cutoff:
        显著性阈值。小麦 KEGG 注释不充分时，推荐用 pvalue。
    """

    # 清理输入基因
    deg_genes = [str(g).strip() for g in gene_list if str(g).strip()]
    deg_genes = list(dict.fromkeys(deg_genes))

    background = load_gene2ko(gene2ko_path)
    pathway_to_kos, pathway_to_name, ko2pathway, pathway2name = load_kegg_mapping(
        ko2pathway_path,
        pathway2name_path
    )

    universe_kos = set(background["KO"].unique())

    deg_background = background[background["gene"].isin(deg_genes)].copy()
    deg_kos = set(deg_background["KO"].unique())

    mapped_gene_count = len(set(deg_genes) & set(background["gene"]))

    summary_df = pd.DataFrame([
        ["输入 DEG 基因数", len(deg_genes)],
        ["成功映射到 KO 的 DEG 基因数", mapped_gene_count],
        ["前景 KO 数", len(deg_kos)],
        ["背景基因数（有 KO 注释）", background["gene"].nunique()],
        ["背景 gene-KO 映射条目数", len(background)],
        ["背景 KO 数", len(universe_kos)],
        ["本地 KEGG pathway 数", len(pathway_to_kos)],
        ["pvalue 阈值", pvalue_cutoff],
        ["最小 pathway KO 数", min_size],
        ["最大 pathway KO 数", max_size],
    ], columns=["项目", "数值"])

    if len(deg_kos) < 3:
        empty = pd.DataFrame()
        return empty, empty, summary_df

    M = len(universe_kos)
    N = len(deg_kos)

    results = []

    for pathway, pathway_kos_all in pathway_to_kos.items():

        pathway_kos_in_bg = set(pathway_kos_all) & universe_kos
        n = len(pathway_kos_in_bg)

        if n < min_size or n > max_size:
            continue

        overlap_kos = pathway_kos_in_bg & deg_kos
        k = len(overlap_kos)

        if k == 0:
            continue

        pvalue = hypergeom.sf(k - 1, M, n, N)

        gene_ratio = f"{k}/{N}"
        bg_ratio = f"{n}/{M}"

        rich_factor = k / n
        fold_enrichment = (k / N) / (n / M)

        results.append({
            "ID": pathway,
            "Description": pathway_to_name.get(pathway, pathway),
            "GeneRatio": gene_ratio,
            "BgRatio": bg_ratio,
            "pvalue": pvalue,
            "geneID": "/".join(sorted(overlap_kos)),
            "Count": k,
            "InputNumber": k,
            "PathwayKoNumberInBackground": n,
            "RichFactor": rich_factor,
            "FoldEnrichment": fold_enrichment
        })

    result_df = pd.DataFrame(results)

    if result_df.empty:
        return result_df, result_df, summary_df

    result_df["p.adjust"] = bh_adjust(result_df["pvalue"].values)
    result_df["qvalue"] = result_df["p.adjust"]

    result_df["minus_log10_pvalue"] = -np.log10(
        result_df["pvalue"].replace(0, 1e-300)
    )

    result_df = result_df.sort_values(
        ["pvalue", "RichFactor", "Count"],
        ascending=[True, False, False]
    )

    result_df = result_df[
        [
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "qvalue",
            "minus_log10_pvalue",
            "geneID",
            "Count",
            "InputNumber",
            "PathwayKoNumberInBackground",
            "RichFactor",
            "FoldEnrichment"
        ]
    ]

    result_df = add_gene_name_columns(
        result_df=result_df,
        background=background,
        deg_genes=deg_genes
    )

    sig_df = result_df[result_df["pvalue"] <= pvalue_cutoff].copy()

    summary_extra = pd.DataFrame([
        ["富集到的 pathway 总数", len(result_df)],
        [f"pvalue <= {pvalue_cutoff} 的 pathway 数", len(sig_df)],
    ], columns=["项目", "数值"])

    summary_df = pd.concat([summary_df, summary_extra], ignore_index=True)

    return result_df, sig_df, summary_df


# ==============================
# 绘图数据准备
# ==============================

def prepare_kegg_plot_df(
    df,
    top_n=20,
    clip_minus_log10_p=True,
    clip_mode="quantile",
    clip_quantile=0.95,
    clip_fixed_value=30,
    min_x_cap=5,
    label_wrap_width=42
):
    """
    准备 KEGG 绘图数据。
    """

    if df.empty:
        return df, 1

    plot_df = df.copy()

    plot_df = plot_df.sort_values(
        ["pvalue", "RichFactor", "Count"],
        ascending=[True, False, False]
    ).head(top_n)

    plot_df["minus_log10_pvalue_raw"] = -np.log10(
        plot_df["pvalue"].replace(0, 1e-300)
    )

    x_cap = get_x_cap(
        plot_df["minus_log10_pvalue_raw"].values,
        clip=clip_minus_log10_p,
        mode=clip_mode,
        quantile=clip_quantile,
        fixed_value=clip_fixed_value,
        min_x_cap=min_x_cap
    )

    if clip_minus_log10_p:
        plot_df["minus_log10_pvalue_plot"] = plot_df["minus_log10_pvalue_raw"].clip(upper=x_cap)
    else:
        plot_df["minus_log10_pvalue_plot"] = plot_df["minus_log10_pvalue_raw"]

    plot_df["Description_wrapped"] = plot_df["Description"].apply(
        lambda x: wrap_text(x, width=label_wrap_width)
    )

    # 反向显示，最显著在上面
    plot_df = plot_df.iloc[::-1].copy()

    return plot_df, x_cap


# ==============================
# KEGG 气泡图
# ==============================

def create_kegg_bubbleplot_bytes(
    df,
    top_n=20,
    bubble_min_size=25,
    bubble_max_size=220,
    clip_minus_log10_p=True,
    clip_mode="quantile",
    clip_quantile=0.95,
    clip_fixed_value=30,
    min_x_cap=5,
    label_wrap_width=42
):
    """
    创建 KEGG 气泡图 PNG bytes。
    """

    if df is None or df.empty:
        return None

    plot_df, x_cap = prepare_kegg_plot_df(
        df=df,
        top_n=top_n,
        clip_minus_log10_p=clip_minus_log10_p,
        clip_mode=clip_mode,
        clip_quantile=clip_quantile,
        clip_fixed_value=clip_fixed_value,
        min_x_cap=min_x_cap,
        label_wrap_width=label_wrap_width
    )

    sizes = scale_bubble_size(
        plot_df["InputNumber"].values,
        min_size=bubble_min_size,
        max_size=bubble_max_size
    )

    height = max(6, 0.45 * len(plot_df))

    fig = plt.figure(figsize=(7.8, height))

    gs = fig.add_gridspec(
        nrows=5,
        ncols=2,
        width_ratios=[1.0, 0.055],
        height_ratios=[0.14, 0.46, 0.06, 0.20, 0.14],
        wspace=0.08,
        hspace=0.05
    )

    ax = fig.add_subplot(gs[:, 0])
    cax = fig.add_subplot(gs[1, 1])
    lax = fig.add_subplot(gs[3, 1])
    lax.axis("off")

    ax.set_facecolor("#EBEBEB")
    fig.patch.set_facecolor("white")

    scatter = ax.scatter(
        plot_df["minus_log10_pvalue_plot"],
        plot_df["Description_wrapped"],
        s=sizes,
        c=plot_df["RichFactor"],
        cmap="coolwarm",
        alpha=0.95,
        edgecolors="black",
        linewidths=0.45
    )

    ax.grid(True, color="white", linewidth=1.2)
    ax.set_axisbelow(True)

    ax.set_xlabel(r"$-\log_{10}(P\mathrm{-value})$", fontsize=12)
    ax.set_ylabel("")

    xmax = max(plot_df["minus_log10_pvalue_plot"].max(), 1)
    ax.set_xlim(0, xmax * 1.10)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(axis="both", labelsize=10)

    cbar = fig.colorbar(scatter, cax=cax)
    cbar.set_label("Rich factor", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    counts = plot_df["InputNumber"].astype(int)

    legend_counts = sorted(set([
        int(counts.min()),
        int(counts.median()),
        int(counts.max())
    ]))

    legend_sizes = scale_bubble_size(
        legend_counts,
        min_size=bubble_min_size,
        max_size=bubble_max_size
    )

    handles = [
        lax.scatter(
            [],
            [],
            s=s,
            facecolors="black",
            edgecolors="black",
            linewidths=0.45,
            alpha=0.95
        )
        for s in legend_sizes
    ]

    lax.legend(
        handles,
        [str(c) for c in legend_counts],
        title="Input\nnumber",
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.82, 0.5),
        scatterpoints=1,
        labelspacing=1.0,
        borderpad=0.2,
        handletextpad=0.8,
        fontsize=8,
        title_fontsize=9
    )

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)

    return buf.getvalue()


# ==============================
# KEGG 条形图
# ==============================

def create_kegg_barplot_bytes(
    df,
    top_n=20,
    clip_minus_log10_p=True,
    clip_mode="quantile",
    clip_quantile=0.95,
    clip_fixed_value=30,
    min_x_cap=5,
    label_wrap_width=42
):
    """
    创建 KEGG 条形图 PNG bytes。
    """

    if df is None or df.empty:
        return None

    plot_df, x_cap = prepare_kegg_plot_df(
        df=df,
        top_n=top_n,
        clip_minus_log10_p=clip_minus_log10_p,
        clip_mode=clip_mode,
        clip_quantile=clip_quantile,
        clip_fixed_value=clip_fixed_value,
        min_x_cap=min_x_cap,
        label_wrap_width=label_wrap_width
    )

    height = max(6, 0.45 * len(plot_df))

    fig, ax = plt.subplots(figsize=(7.2, height))

    ax.set_facecolor("#EBEBEB")
    fig.patch.set_facecolor("white")

    values = plot_df["RichFactor"].values

    if values.max() == values.min():
        norm = plt.Normalize(vmin=0, vmax=max(values.max(), 1e-6))
    else:
        norm = plt.Normalize(vmin=values.min(), vmax=values.max())

    cmap = plt.cm.coolwarm

    ax.barh(
        plot_df["Description_wrapped"],
        plot_df["minus_log10_pvalue_plot"],
        color=cmap(norm(values)),
        edgecolor="none"
    )

    # 竖向白色网格线
    ax.grid(
        True,
        axis="x",
        color="white",
        linewidth=1.2,
        alpha=0.95
    )

    # 横向白色纹理线，让条形图更接近气泡图风格
    ax.grid(
        True,
        axis="y",
        color="white",
        linewidth=0.8,
        alpha=0.85
    )

    ax.set_axisbelow(True)

    ax.set_xlabel(r"$-\log_{10}(P\mathrm{-value})$", fontsize=12)
    ax.set_ylabel("")

    xmax = max(plot_df["minus_log10_pvalue_plot"].max(), 1)
    ax.set_xlim(0, xmax * 1.12)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(axis="both", labelsize=10)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cbar = plt.colorbar(
        sm,
        ax=ax,
        fraction=0.08,
        pad=0.035,
        shrink=0.5,
        aspect=13.85
    )
    cbar.set_label("Rich factor", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)

    return buf.getvalue()
