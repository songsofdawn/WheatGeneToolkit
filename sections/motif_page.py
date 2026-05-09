from pathlib import Path

import pandas as pd
import streamlit as st

from utils.jaspar_pwm_scan import (
    load_jaspar_pwm,
    parse_fasta_or_plain,
    scan_sequences_with_jaspar_significance,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
JASPAR_PWM_PATH = PROJECT_ROOT / "data" / "motif_db" / "jaspar_plants" / "jaspar_plants_pwm.json"
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
    "p_value",
    "q_value",
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
    "total_significant_hits",
    "best_q_value",
    "max_relative_score",
    "best_matched_seq",
    "best_position",
]


@st.cache_data(show_spinner=False)
def _load_cached_motifs(json_path: str):
    """
    缓存 JASPAR PWM JSON，避免 Streamlit 每次交互都重新读取文件。
    """
    return load_jaspar_pwm(json_path)


def _filter_motifs(motifs, keyword: str):
    """
    按 matrix_id 或 TF name 关键词筛选 motif；为空时返回全部。
    """
    keyword = (keyword or "").strip().lower()
    if not keyword:
        return motifs

    filtered = []
    for motif in motifs:
        matrix_id = str(motif.get("matrix_id", "")).lower()
        tf_name = str(motif.get("name", "") or (motif.get("metadata") or {}).get("name", "")).lower()
        if keyword in matrix_id or keyword in tf_name:
            filtered.append(motif)
    return filtered


def _build_significant_summary(significant_df: pd.DataFrame) -> pd.DataFrame:
    """
    对显著 TFBS 结果按序列和 motif 分组汇总。
    """
    if significant_df.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    group_cols = ["sequence_id", "matrix_id", "tf_name", "consensus"]
    summary_df = (
        significant_df.groupby(group_cols, dropna=False)
        .agg(
            total_significant_hits=("sequence_id", "size"),
            best_q_value=("q_value", "min"),
            max_relative_score=("relative_score", "max"),
        )
        .reset_index()
    )

    best_rows = (
        significant_df.sort_values(["q_value", "relative_score"], ascending=[True, False])
        .drop_duplicates(group_cols)
        .copy()
    )
    best_rows["best_position"] = (
        best_rows["start"].astype(str) + "-" + best_rows["end"].astype(str) + "(" + best_rows["strand"].astype(str) + ")"
    )
    best_rows = best_rows[group_cols + ["matched_seq", "best_position"]].rename(
        columns={"matched_seq": "best_matched_seq"}
    )

    summary_df = summary_df.merge(best_rows, on=group_cols, how="left")
    return summary_df[SUMMARY_COLUMNS].sort_values(
        ["best_q_value", "total_significant_hits", "max_relative_score"],
        ascending=[True, False, False],
    )


def render():
    st.header("启动子 JASPAR Plants PWM motif 分析")

    st.markdown(
        """
        本模块使用 JASPAR CORE Plants non-redundant PFMs 转换得到的 PWM，
        用于预测输入启动子中的潜在转录因子结合位点。

        结果是序列层面的预测，不等同于真实结合或真实调控证据。
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

    st.info(f"已加载 JASPAR Plants PWM motifs: {len(motifs)}")

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
            value=0.85,
            step=0.01,
        )
    with col2:
        qvalue_cutoff = st.number_input(
            "q-value cutoff",
            min_value=0.001,
            max_value=1.0,
            value=0.05,
            step=0.01,
            format="%.3f",
        )
    with col3:
        scan_reverse = st.checkbox("同时扫描反向互补链", value=True)

    col4, col5, col6 = st.columns(3)
    with col4:
        max_total_hits = st.number_input(
            "最大输出 hits 数",
            min_value=100,
            max_value=100000,
            value=20000,
            step=1000,
        )
    with col5:
        background_mode = st.selectbox(
            "背景模型",
            options=["input", "uniform"],
            format_func=lambda value: {
                "input": "使用输入序列估计背景",
                "uniform": "均匀背景 A/C/G/T=0.25",
            }[value],
        )
    with col6:
        n_background_samples = st.selectbox(
            "背景模拟次数",
            options=[10000, 20000, 50000, 100000],
            index=2,
        )

    st.caption("背景模拟次数越高，p-value 估计越稳定，但扫描会更慢；测试时可先选择 10000。")

    with st.expander("高级筛选", expanded=False):
        motif_keyword = st.text_input(
            "按 matrix_id 或 TF name 关键词筛选 motif",
            placeholder="例如 MA0561 或 PIF4；留空表示扫描全部 motif",
        )

    selected_motifs = _filter_motifs(motifs, motif_keyword)
    st.caption(f"当前参与扫描的 motifs: {len(selected_motifs)}")

    if st.button("开始 JASPAR PWM 扫描", key="btn_jaspar_pwm_scan"):
        records = parse_fasta_or_plain(sequence_text)
        records = {name: seq for name, seq in records.items() if seq}

        if not records:
            st.warning("请粘贴 FASTA 或 DNA 序列后再开始扫描。")
            st.stop()

        if not selected_motifs:
            st.warning("当前筛选条件下没有可扫描的 motif，请调整关键词。")
            st.stop()

        total_length = sum(len(seq) for seq in records.values())
        st.info(f"输入序列数量: {len(records)}；总长度: {total_length} bp")
        if total_length > 100000:
            st.warning("输入序列总长度超过 100000 bp，PWM 全库扫描和背景模拟可能较慢。可提高初筛阈值、降低模拟次数或使用高级筛选减少 motif 数量。")

        with st.spinner("正在进行 JASPAR PWM 扫描和背景显著性估计，请稍候..."):
            result_df = scan_sequences_with_jaspar_significance(
                records=records,
                motifs=selected_motifs,
                relative_cutoff=relative_cutoff,
                qvalue_cutoff=qvalue_cutoff,
                scan_reverse=scan_reverse,
                selected_matrix_ids=None,
                max_total_hits=int(max_total_hits),
                background_mode=background_mode,
                n_background_samples=int(n_background_samples),
                random_seed=123,
            )

        if result_df.empty:
            st.warning("没有发现超过 relative score 初筛阈值的潜在 TF binding sites。可以适当降低初筛阈值，例如 0.80。")
            st.stop()

        result_df = result_df[DETAIL_COLUMNS].sort_values(["q_value", "relative_score"], ascending=[True, False])

        if len(result_df) >= int(max_total_hits):
            st.warning("候选结果已达到最大输出 hits 数限制，实际命中数可能更多。可提高初筛阈值或增加最大输出 hits 数。")

        significant_df = result_df[result_df["significant"] == True].copy()
        significant_df = significant_df.sort_values(["q_value", "relative_score"], ascending=[True, False])

        st.success(
            f"扫描完成：relative score 初筛候选 {len(result_df)} 个；"
            f"q-value <= {qvalue_cutoff:.3f} 的显著候选 {len(significant_df)} 个。"
        )

        st.subheader("显著候选 TFBS")
        if significant_df.empty:
            st.warning(
                "没有发现 q-value 小于阈值的显著候选 TFBS。可以尝试降低 relative score 初筛阈值，"
                "或增加背景模拟次数，但请谨慎解释。"
            )
        else:
            st.dataframe(significant_df, use_container_width=True)
            st.download_button(
                "下载显著结果 CSV",
                data=significant_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_significant_hits.csv",
                mime="text/csv",
            )

            summary_df = _build_significant_summary(significant_df)
            st.subheader("显著 motif 命中汇总")
            st.dataframe(summary_df, use_container_width=True)
            st.download_button(
                "下载显著汇总结果 CSV",
                data=summary_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_significant_summary.csv",
                mime="text/csv",
            )

        with st.expander("查看所有 relative score 初筛候选 hit", expanded=False):
            st.dataframe(result_df, use_container_width=True)
            st.download_button(
                "下载全部候选结果 CSV",
                data=result_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_candidate_hits.csv",
                mime="text/csv",
            )

    st.warning(
        "注意：relative_score 表示窗口与 PWM 最佳模式的相似度；p-value 表示在背景模型下获得不低于当前 PWM score 的概率；"
        "q-value 是对多个 motif 和多个窗口检验后的 FDR 校正结果。默认主结果只报告 q-value <= 0.05 的候选 TFBS。"
        "即使 q-value 显著，也仍然是序列预测，不等同于真实结合证据。"
    )

    st.info(
        "本模块采用 PWM score 对启动子窗口进行打分，并通过随机背景模型估计 p-value，再使用 Benjamini-Hochberg "
        "方法进行多重检验校正得到 q-value。主结果默认仅展示 q-value 小于阈值的候选 TFBS。"
        "该策略相比单纯 relative score 阈值可以减少短核心 motif 的大量假阳性，但结果仍需结合表达数据、ATAC-seq、"
        "ChIP-seq 或实验验证。"
    )
