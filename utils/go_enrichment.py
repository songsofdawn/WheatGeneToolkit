import io
import math
from functools import lru_cache
from pathlib import Path
from typing import Set, Dict, List, Tuple, Optional

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
    top_n: int = 10,
    plot_metric: str = "fdr",
    figsize=(12, 12)
) -> Optional[bytes]:
    if df.empty:
        return None

    df = df[df["ontology"].isin(["BP", "CC", "MF"])].copy()
    if df.empty:
        return None

    if plot_metric == "fdr":
        score_col = "minus_log10_padj"
        x_label = "-log10(FDR)"
        sort_col = "p.adjust"
    elif plot_metric == "pvalue":
        score_col = "minus_log10_pvalue"
        x_label = "-log10(P-value)"
        sort_col = "pvalue"
    else:
        raise ValueError("plot_metric 只能是 'fdr' 或 'pvalue'")

    subsets = []
    for onto in ["BP", "CC", "MF"]:
        sub = df[df["ontology"] == onto].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(sort_col, ascending=True).head(top_n)
        sub["Description_short"] = sub["go_term_name"].fillna(sub["go_id"])
        subsets.append((onto, sub))

    if not subsets:
        return None

    n_panels = len(subsets)
    fig, axes = plt.subplots(n_panels, 1, figsize=figsize, constrained_layout=True)

    if n_panels == 1:
        axes = [axes]

    for ax, (onto, sub) in zip(axes, subsets):
        sub = sub.sort_values(score_col, ascending=True)
        ax.barh(sub["Description_short"], sub[score_col])
        ax.set_title(f"{onto} enrichment")
        ax.set_xlabel(x_label)
        ax.set_ylabel("GO term")

    buf = io.BytesIO()
    plt.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.getvalue()
