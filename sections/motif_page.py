from pathlib import Path

import streamlit as st

from utils.jaspar_pwm_scan import (
    load_jaspar_pwm,
    parse_fasta_or_plain,
    scan_sequences_with_jaspar,
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
    "distance_to_sequence_end",
    "species",
    "tax_group",
    "family",
    "class",
    "collection",
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


def _build_summary(result_df):
    """
    生成 motif 命中汇总表。
    """
    return (
        result_df.groupby(["matrix_id", "tf_name", "consensus"], dropna=False)
        .agg(
            sequence_count=("sequence_id", "nunique"),
            total_hits=("sequence_id", "size"),
            max_relative_score=("relative_score", "max"),
            mean_relative_score=("relative_score", "mean"),
        )
        .reset_index()
        .sort_values(["total_hits", "max_relative_score"], ascending=[False, False])
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
        cutoff = st.slider(
            "relative score cutoff",
            min_value=0.70,
            max_value=0.99,
            value=0.85,
            step=0.01,
        )
    with col2:
        scan_reverse = st.checkbox("同时扫描反向互补链", value=True)
    with col3:
        max_total_hits = st.number_input(
            "最大输出 hits 数",
            min_value=100,
            max_value=100000,
            value=20000,
            step=1000,
        )

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
            st.warning("输入序列总长度超过 100000 bp，PWM 全库扫描可能较慢。可提高 cutoff 或使用高级筛选减少 motif 数量。")

        with st.spinner("正在进行 JASPAR PWM 扫描，请稍候..."):
            result_df = scan_sequences_with_jaspar(
                records=records,
                motifs=selected_motifs,
                cutoff=cutoff,
                scan_reverse=scan_reverse,
                selected_matrix_ids=None,
                max_total_hits=int(max_total_hits),
            )

        if result_df.empty:
            st.warning("没有发现超过阈值的潜在 TF binding sites。可以适当降低 cutoff，例如 0.80。")
        else:
            # 固定列顺序，便于页面查看和 CSV 后续分析。
            result_df = result_df[DETAIL_COLUMNS]
            summary_df = _build_summary(result_df)

            if len(result_df) >= int(max_total_hits):
                st.warning("结果已达到最大输出 hits 数限制，实际命中数可能更多。可提高 cutoff 或增加最大输出 hits 数。")

            st.success(f"扫描完成，共发现 {len(result_df)} 个潜在 TF binding sites。")

            st.subheader("潜在 TF binding sites 明细")
            st.dataframe(result_df, use_container_width=True)
            st.download_button(
                "下载明细结果 CSV",
                data=result_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_hits.csv",
                mime="text/csv",
            )

            st.subheader("motif 命中汇总")
            st.dataframe(summary_df, use_container_width=True)
            st.download_button(
                "下载汇总结果 CSV",
                data=summary_df.to_csv(index=False).encode("utf-8-sig"),
                file_name="jaspar_plants_pwm_summary.csv",
                mime="text/csv",
            )

    st.warning(
        "注意：PWM 命中只表示启动子中存在与该转录因子结合偏好相似的序列片段，"
        "不等同于真实调控关系。建议结合表达数据、保守性、ATAC-seq、ChIP-seq "
        "或实验验证进一步确认。"
    )
