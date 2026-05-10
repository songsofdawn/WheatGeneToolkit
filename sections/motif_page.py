from pathlib import Path

import pandas as pd
import streamlit as st

from utils.db_query import get_promoter
from utils.jaspar_pwm_scan import (
    load_jaspar_pwm,
    load_precomputed_thresholds,
    parse_fasta_or_plain,
    p_level_passes_cutoff,
    scan_sequences_with_jaspar_thresholds,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
JASPAR_DIR = PROJECT_ROOT / "data" / "motif_db" / "jaspar_plants"
JASPAR_PWM_PATH = JASPAR_DIR / "jaspar_plants_pwm.json"
THRESHOLD_PATHS = {
    "uniform": JASPAR_DIR / "jaspar_background_thresholds_uniform.json",
    "cs_promoter": JASPAR_DIR / "jaspar_background_thresholds_cs_promoter.json",
    "fielder_promoter": JASPAR_DIR / "jaspar_background_thresholds_fielder_promoter.json",
}
EXAMPLE_TAGW2_GENE_ID = "TraesCS6A02G189300"
EXAMPLE_TAGW2_FASTA_HEADER = ">TaGW2-A|TraesCS6A02G189300|promoter_2000"

MAIN_COLUMNS = [
    "sequence_id",
    "matrix_id",
    "tf_name",
    "consensus",
    "motif_length",
    "start",
    "end",
    "strand",
    "matched_seq",
    "raw_score",
    "relative_score",
    "p_level",
    "confidence_level",
    "species",
    "family",
    "class",
]
DETAIL_COLUMNS = [
    "sequence_id",
    "matrix_id",
    "tf_name",
    "consensus",
    "motif_length",
    "start",
    "end",
    "strand",
    "matched_seq",
    "raw_score",
    "relative_score",
    "p_level",
    "significance_rank",
    "confidence_level",
    "significant",
    "distance_to_sequence_end",
    "species",
    "tax_group",
    "family",
    "class",
    "collection",
]
SUMMARY_COLUMNS = [
    "sequence_id",
    "matrix_id",
    "tf_name",
    "consensus",
    "total_hits",
    "best_p_level",
    "best_confidence_level",
    "max_relative_score",
    "mean_relative_score",
    "best_matched_seq",
    "best_position",
]


@st.cache_data(show_spinner=False)
def _load_cached_motifs(json_path: str):
    """
    缓存 JASPAR PWM JSON，避免 Streamlit 每次交互都重新读取文件。
    """
    return load_jaspar_pwm(json_path)


@st.cache_data(show_spinner=False)
def _load_cached_thresholds(json_path: str):
    """
    缓存轻量级背景阈值表。文件不存在时返回 None。
    """
    return load_precomputed_thresholds(json_path)


def _load_threshold_payloads():
    """
    尝试加载所有 threshold JSON。
    """
    return {key: _load_cached_thresholds(str(path)) for key, path in THRESHOLD_PATHS.items()}


def _wrap_sequence(sequence: str, width: int = 80) -> str:
    """
    将 DNA 序列按固定宽度换行，便于在文本框中以 FASTA 格式展示。
    """
    sequence = "".join(str(sequence or "").split()).upper()
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def _load_tagw2_example_fasta():
    """
    从项目已有 Chinese Spring promoter 查询函数中读取 TaGW2-A 2000 bp 启动子。
    不在代码中硬编码序列，保证示例和本地数据库保持一致。
    """
    try:
        promoter_df = get_promoter(EXAMPLE_TAGW2_GENE_ID)
    except Exception:
        return None

    if promoter_df is None or promoter_df.empty or "promoter_sequence" not in promoter_df.columns:
        return None

    promoter_seq = _wrap_sequence(promoter_df.iloc[0].get("promoter_sequence", ""))
    if not promoter_seq:
        return None

    return f"{EXAMPLE_TAGW2_FASTA_HEADER}\n{promoter_seq}\n"


def _has_total_ic(motifs) -> bool:
    return any("total_ic" in motif for motif in motifs)


def _filter_motifs(motifs, keyword: str, min_motif_length: int = 6, min_total_ic=None):
    """
    按关键词、motif 长度和 total_ic 过滤 motif。
    """
    keyword = (keyword or "").strip().lower()
    filtered = []

    for motif in motifs:
        matrix_id = str(motif.get("matrix_id", "")).lower()
        tf_name = str(motif.get("name", "") or (motif.get("metadata") or {}).get("name", "")).lower()
        motif_length = int(motif.get("length") or 0)

        if keyword and keyword not in matrix_id and keyword not in tf_name:
            continue
        if motif_length < int(min_motif_length):
            continue

        if min_total_ic is not None and "total_ic" in motif:
            try:
                if float(motif.get("total_ic", 0.0)) < float(min_total_ic):
                    continue
            except (TypeError, ValueError):
                continue

        filtered.append(motif)

    return filtered


def _threshold_option_label(option: str, payloads: dict) -> str:
    labels = {
        "uniform": "预计算均匀背景 A/C/G/T=0.25（推荐，速度快）",
        "cs_promoter": "Chinese Spring 启动子经验背景",
        "fielder_promoter": "Fielder 启动子经验背景",
    }
    status = "已加载" if payloads.get(option) else "未找到"
    return f"{labels[option]} - {status}"


def _build_summary(main_df: pd.DataFrame) -> pd.DataFrame:
    """
    对主结果按序列和 motif 分组汇总。
    """
    if main_df.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    group_cols = ["sequence_id", "matrix_id", "tf_name", "consensus"]
    summary_df = (
        main_df.groupby(group_cols, dropna=False)
        .agg(
            total_hits=("sequence_id", "size"),
            best_rank=("significance_rank", "min"),
            max_relative_score=("relative_score", "max"),
            mean_relative_score=("relative_score", "mean"),
        )
        .reset_index()
    )

    best_rows = (
        main_df.sort_values(["significance_rank", "relative_score", "raw_score"], ascending=[True, False, False])
        .drop_duplicates(group_cols)
        .copy()
    )
    best_rows["best_position"] = (
        best_rows["start"].astype(str)
        + "-"
        + best_rows["end"].astype(str)
        + "("
        + best_rows["strand"].astype(str)
        + ")"
    )
    best_rows = best_rows[
        group_cols + ["p_level", "confidence_level", "matched_seq", "best_position"]
    ].rename(
        columns={
            "p_level": "best_p_level",
            "confidence_level": "best_confidence_level",
            "matched_seq": "best_matched_seq",
        }
    )

    summary_df = summary_df.merge(best_rows, on=group_cols, how="left")
    summary_df = summary_df.sort_values(["best_rank", "max_relative_score"], ascending=[True, False])
    return summary_df[SUMMARY_COLUMNS]


def render():
    st.header("启动子 JASPAR Plants PWM motif 分析")

    st.markdown(
        """
        本模块使用 JASPAR Plants PWM 对启动子滑动窗口进行打分，并使用离线预计算的轻量级背景
        score cutoff 表对命中结果进行显著性分级。该模式不加载完整背景分布 JSON，适合 Streamlit 快速运行。
        """
    )

    if not JASPAR_PWM_PATH.exists():
        st.error("未找到 JASPAR Plants PWM JSON 文件，请检查 data/motif_db/jaspar_plants/ 目录。")
        st.code(str(JASPAR_PWM_PATH))
        st.stop()

    try:
        motifs = _load_cached_motifs(str(JASPAR_PWM_PATH))
    except Exception as exc:
        st.error("读取 JASPAR Plants PWM JSON 失败。")
        st.exception(exc)
        st.stop()

    threshold_payloads = _load_threshold_payloads()
    loaded_thresholds = [key for key, payload in threshold_payloads.items() if payload]
    if loaded_thresholds:
        st.success("已加载轻量级背景阈值表，适合 Streamlit 快速分析。")
    else:
        st.warning(
            "未找到轻量级背景阈值表。请先运行：python scripts/build_jaspar_background.py。"
        )

    st.info(f"已加载 JASPAR Plants PWM motifs: {len(motifs)}")

    st.subheader("加载示例")
    st.caption(
        "TaGW2-A 是小麦籽粒大小和粒重相关基因。该示例用于展示 JASPAR Plants PWM motif "
        "分析流程，结果仅代表潜在 TFBS 预测。"
    )
    if st.button("加载示例：TaGW2-A / TraesCS6A02G189300", key="btn_load_tagw2_motif_example"):
        example_fasta = _load_tagw2_example_fasta()
        if example_fasta is None:
            st.error("未能从本地启动子数据库中找到 TraesCS6A02G189300，请确认 Chinese Spring promoter 数据库是否存在。")
        else:
            st.session_state["jaspar_promoter_input"] = example_fasta
            st.success("已加载 TaGW2-A / TraesCS6A02G189300 的 2000 bp 启动子示例。")

    sequence_text = st.text_area(
        "粘贴启动子 FASTA 或纯 DNA 序列",
        height=220,
        placeholder=">test_promoter\nATGCGTACGTACGTACGTGACGTAGCTAGCTA",
        key="jaspar_promoter_input",
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        relative_cutoff = st.slider(
            "relative score 初筛阈值",
            min_value=0.70,
            max_value=0.99,
            value=0.90,
            step=0.01,
        )
    with col2:
        p_level_cutoff = st.selectbox(
            "p-level 阈值",
            options=[0.05, 0.01, 0.001, 0.0001],
            index=1,
            format_func=lambda value: f"p_level <= {value:g}",
        )
    with col3:
        scan_reverse = st.checkbox("同时扫描反向互补链", value=True)

    col4, col5 = st.columns(2)
    with col4:
        max_total_hits = st.number_input(
            "最大候选 hits 数",
            min_value=100,
            max_value=100000,
            value=20000,
            step=1000,
        )
    with col5:
        background_choice = st.selectbox(
            "背景阈值表",
            options=["uniform", "cs_promoter", "fielder_promoter"],
            index=0,
            format_func=lambda option: _threshold_option_label(option, threshold_payloads),
        )

    selected_thresholds = threshold_payloads.get(background_choice)
    if selected_thresholds:
        st.caption(
            f"当前阈值表样本数: {selected_thresholds.get('n_samples', '未知')}；"
            f"p-levels: {selected_thresholds.get('p_levels', [])}"
        )
    else:
        st.error("当前选择的阈值表不存在。请先生成 threshold JSON 后再扫描。")

    has_total_ic = _has_total_ic(motifs)
    with st.expander("高级筛选", expanded=False):
        motif_keyword = st.text_input(
            "按 matrix_id 或 TF name 关键词筛选 motif",
            placeholder="例如 MA0561 或 PIF4；留空表示扫描全部 motif",
        )
        col_a, col_b, col_c = st.columns(3)
        with col_a:
            min_motif_length = st.number_input(
                "最小 motif_length",
                min_value=1,
                max_value=50,
                value=6,
                step=1,
            )
        with col_b:
            if has_total_ic:
                min_total_ic = st.number_input(
                    "最小 total_ic",
                    min_value=0.0,
                    max_value=50.0,
                    value=6.0,
                    step=0.5,
                )
            else:
                min_total_ic = None
                st.caption("当前 PWM JSON 不含 total_ic 字段，仅使用 motif_length 过滤。")
        with col_c:
            top_n_per_motif_sequence = st.number_input(
                "每个 motif 每条序列 top hits",
                min_value=1,
                max_value=20,
                value=3,
                step=1,
            )

    selected_motifs = _filter_motifs(
        motifs,
        keyword=motif_keyword,
        min_motif_length=int(min_motif_length),
        min_total_ic=min_total_ic,
    )
    st.caption(f"当前参与扫描的 motifs: {len(selected_motifs)}")

    if st.button("开始 JASPAR PWM 扫描", key="btn_jaspar_pwm_scan"):
        if not selected_thresholds:
            st.error("没有可用的轻量级背景阈值表。请先运行 python scripts/build_jaspar_background.py。")
            st.stop()

        records = parse_fasta_or_plain(sequence_text)
        records = {name: seq for name, seq in records.items() if seq}

        if not records:
            st.warning("请粘贴 FASTA 或 DNA 序列后再开始扫描。")
            st.stop()

        if not selected_motifs:
            st.warning("当前筛选条件下没有可扫描的 motif，请调整关键词、motif_length 或 total_ic。")
            st.stop()

        total_length = sum(len(seq) for seq in records.values())
        st.info(f"输入序列数量: {len(records)}；总长度: {total_length} bp")
        if total_length > 100000:
            st.warning("输入序列总长度超过 100000 bp，PWM 扫描可能较慢。可提高初筛阈值或使用高级筛选减少 motif 数量。")

        with st.spinner("正在进行 JASPAR PWM 扫描和 p-level 显著性分级，请稍候..."):
            candidate_df = scan_sequences_with_jaspar_thresholds(
                records=records,
                motifs=selected_motifs,
                relative_cutoff=relative_cutoff,
                p_level_cutoff=p_level_cutoff,
                scan_reverse=scan_reverse,
                selected_matrix_ids=None,
                max_total_hits=int(max_total_hits),
                precomputed_thresholds=selected_thresholds,
                top_n_per_motif_sequence=int(top_n_per_motif_sequence),
            )

        if candidate_df.empty:
            st.warning("没有发现超过 relative score 初筛阈值的潜在 TF binding sites。可以适当降低初筛阈值，例如 0.85。")
            st.stop()

        candidate_df = candidate_df[DETAIL_COLUMNS].sort_values(
            ["significance_rank", "relative_score", "raw_score"],
            ascending=[True, False, False],
        )
        main_df = candidate_df[
            candidate_df["p_level"].apply(lambda value: p_level_passes_cutoff(value, p_level_cutoff))
        ].copy()

        st.success(
            f"扫描完成：候选 hits {len(candidate_df)} 个；"
            f"通过 p_level <= {p_level_cutoff:g} 的主结果 {len(main_df)} 个。"
        )

        st.subheader("显著主结果表")
        if main_df.empty:
            st.warning(
                "没有发现满足当前 p-level 阈值的候选 TFBS。可以尝试降低 relative score 初筛阈值，"
                "或使用更宽松的 p-level 阈值，但请谨慎解释。"
            )
        else:
            st.dataframe(main_df[MAIN_COLUMNS], use_container_width=True)
            st.download_button(
                "下载显著结果 CSV",
                data=main_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_threshold_hits.csv",
                mime="text/csv",
            )

            summary_df = _build_summary(main_df)
            st.subheader("汇总表")
            st.dataframe(summary_df, use_container_width=True)
            st.download_button(
                "下载汇总结果 CSV",
                data=summary_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_threshold_summary.csv",
                mime="text/csv",
            )

        with st.expander("候选结果表：查看所有 relative score 初筛通过的 top hits", expanded=False):
            st.dataframe(candidate_df, use_container_width=True)
            st.download_button(
                "下载全部候选结果 CSV",
                data=candidate_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_threshold_candidates.csv",
                mime="text/csv",
            )

    st.info(
        "本模块使用预计算背景 score cutoff 对 PWM 命中进行显著性分级。p_level 表示该 raw_score 至少达到某个背景显著性水平。"
        "例如 p_level <= 0.001 表示随机背景中约 0.1% 的窗口能达到或超过该分数。"
        "该方法比单纯 relative_score 更能减少短核心 motif 假阳性，同时保持 Streamlit 快速运行。"
        "由于该模式是分级近似，不输出精确 p-value/q-value。如需精确统计检验，建议使用离线 FIMO/MEME Suite 或后续高级离线模式。"
        "即使分级显著，结果仍然只是序列层面的结合位点预测，不等同于真实调控关系；建议结合表达数据、ATAC-seq、ChIP-seq、"
        "保守性分析或实验验证进一步确认。"
    )
