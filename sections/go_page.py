import io
import os
import time

import pandas as pd
import streamlit as st

from app_shared import (
    get_go_example_text,
    get_go_mapping_paths,
    preprocess_gene_ids,
    render_example_tools,
    run_cached_go_enrichment,
    show_dataframe_preview,
    show_input_cleanup_notice,
    show_large_input_notice,
)
from utils.go_enrichment import create_go_barplot_panel_bytes, create_go_bubbleplot_bytes


def _df_cache_text(df):
    return df.to_csv(index=False)


@st.cache_data(show_spinner=False)
def _cached_go_barplot_bytes(df_text, ontology, top_n, qvalue_col, count_col, label_wrap_width):
    df = pd.read_csv(io.StringIO(df_text))
    return create_go_barplot_panel_bytes(
        df=df,
        ontology=ontology,
        top_n=top_n,
        qvalue_col=qvalue_col,
        count_col=count_col,
        label_wrap_width=label_wrap_width,
    )


@st.cache_data(show_spinner=False)
def _cached_go_bubbleplot_bytes(df_text, ontology, top_n, qvalue_col, count_col, label_wrap_width):
    df = pd.read_csv(io.StringIO(df_text))
    return create_go_bubbleplot_bytes(
        df=df,
        ontology=ontology,
        top_n=top_n,
        qvalue_col=qvalue_col,
        count_col=count_col,
        label_wrap_width=label_wrap_width,
    )


def render():
    st.header("GO 富集分析")
    st.caption("输入 DEG 列表（一行一个基因号），输出 GO 富集结果表、条形图和气泡图。")

    uploaded_file = st.file_uploader("上传 DEG TXT 文件（一行一个基因号）", type=["txt"], key="file_go_deg")
    render_example_tools(
        input_key="input_go_deg",
        example_text=get_go_example_text(),
        load_label="加载 GO 示例 DEG",
        download_label="下载 GO 示例 DEG 文件",
        file_name="example_go_deg_genes.txt",
        help_text="示例 DEG 从 data/example_go_deg_genes.txt 读取；你可以直接修改该文件来更新 GO 示例。",
    )
    manual_input = st.text_area("或者手动输入 DEG（每行一个）", height=200, key="input_go_deg")

    gene_ids, cleanup_info = preprocess_gene_ids(uploaded_file, manual_input)
    if not gene_ids:
        st.info("请上传 DEG 文件或手动输入基因号")
        st.stop()
    show_input_cleanup_notice(cleanup_info)
    show_large_input_notice(len(gene_ids), task_name="GO 富集分析", threshold=2000)

    col1, col2, col3 = st.columns(3)
    with col1:
        top_n = st.number_input("每个大类显示前 N 个 term", min_value=5, max_value=30, value=15, step=1)
    with col2:
        padj_cutoff = st.number_input("FDR 阈值", min_value=0.0001, max_value=1.0, value=0.05, step=0.01, format="%.4f")
    with col3:
        plot_type = st.radio(
            "图形类型",
            ["Bar plot", "Bubble plot", "Both"],
            index=0,
            horizontal=True,
        )

    with st.expander("绘图参数"):
        label_wrap_width = st.number_input("GO term 换行宽度", min_value=20, max_value=70, value=35, step=5)
        st.caption("GO 富集图使用 -log₁₀(q-value) 映射颜色；数值越大表示 q-value 越小，颜色越偏红。颜色范围会根据当前 BP/CC/MF 图中展示的 term 自适应调整，因此不同图之间的色带范围可能略有不同，主要用于比较同一张图内部的相对显著性。")

    min_size = st.number_input("最小 GO 基因集大小", min_value=1, max_value=50, value=3, step=1)
    max_size = st.number_input("最大 GO 基因集大小", min_value=10, max_value=10000, value=2000, step=10)

    go_paths = get_go_mapping_paths()
    term2gene_path = go_paths["term2gene"]
    term2name_path = go_paths["term2name"]
    metadata_path = go_paths["metadata"]
    background_path = go_paths["background"]

    required_files = [term2gene_path, term2name_path, metadata_path, background_path]
    missing_files = [fp for fp in required_files if not os.path.exists(fp)]
    if missing_files:
        st.error("以下 GO 注释文件不存在，请检查 data/go_mapping/ 目录：")
        for fp in missing_files:
            st.code(fp)
        st.stop()

    st.info(f"待分析基因数: {len(gene_ids)}")

    if st.button("开始 GO 富集分析", key="btn_go_enrichment"):
        started_at = time.perf_counter()
        with st.spinner("正在进行 GO 富集分析，请稍候..."):
            try:
                results_df, sig_df, summary_df = run_cached_go_enrichment(
                    gene_ids=tuple(gene_ids),
                    term2gene_path=term2gene_path,
                    term2name_path=term2name_path,
                    metadata_path=metadata_path,
                    background_path=background_path,
                    min_size=min_size,
                    max_size=max_size,
                    padj_cutoff=padj_cutoff,
                )

                if results_df.empty:
                    st.warning("没有得到任何 GO 富集结果，请检查输入基因 ID 是否与背景一致。")
                    st.stop()

                st.success("GO 富集分析完成")
                st.caption(f"富集分析用时：{time.perf_counter() - started_at:.2f} 秒")
                st.subheader("分析摘要")
                show_dataframe_preview(summary_df, label="GO 分析摘要", key="show_all_go_summary")

                st.subheader("显著富集结果")
                if sig_df.empty:
                    st.warning("当前 FDR 阈值下没有显著富集条目，下面展示全部结果。")
                    show_dataframe_preview(results_df, label="GO 全部富集结果", key="show_all_go_results_no_sig")
                else:
                    show_dataframe_preview(sig_df, label="GO 显著富集结果", key="show_all_go_sig")

                plot_df = sig_df.copy()
                if plot_df.empty:
                    plot_df = results_df.copy()
                plot_df_text = _df_cache_text(plot_df)

                for ontology, label in [("BP", "BP"), ("CC", "CC"), ("MF", "MF")]:
                    if plot_df.empty or plot_df[plot_df["ontology"] == ontology].empty:
                        st.info(f"{label} 类别未检测到显著富集的 GO term。")

                show_bar = plot_type in {"Bar plot", "Both"}
                show_bubble = plot_type in {"Bubble plot", "Both"}

                if show_bar:
                    plot_started_at = time.perf_counter()
                    st.subheader("GO 富集条形图")
                    for ontology in ["BP", "CC", "MF"]:
                        barplot_bytes = _cached_go_barplot_bytes(
                            plot_df_text,
                            ontology,
                            top_n,
                            "qvalue",
                            "Count",
                            label_wrap_width,
                        )
                        if barplot_bytes is None:
                            st.info(f"{ontology} 类别没有可绘制的 GO 富集条形图。")
                            continue

                        st.image(
                            barplot_bytes,
                            caption=f"{ontology} GO enrichment barplot",
                            use_container_width=True,
                        )
                        st.download_button(
                            f"下载 {ontology} GO 条形图 PNG",
                            data=barplot_bytes,
                            file_name=f"GO_{ontology}_enrichment_barplot.png",
                            mime="image/png",
                            key=f"download_go_{ontology}_barplot",
                        )
                    st.caption(f"GO 条形图生成用时：{time.perf_counter() - plot_started_at:.2f} 秒")

                if show_bubble:
                    plot_started_at = time.perf_counter()
                    st.subheader("GO 富集气泡图")
                    for ontology in ["BP", "CC", "MF"]:
                        bubbleplot_bytes = _cached_go_bubbleplot_bytes(
                            plot_df_text,
                            ontology,
                            top_n,
                            "qvalue",
                            "Count",
                            label_wrap_width,
                        )
                        if bubbleplot_bytes is None:
                            st.info(f"{ontology} 类别没有可绘制的 GO 富集气泡图。")
                            continue

                        st.image(
                            bubbleplot_bytes,
                            caption=f"{ontology} GO bubble plot",
                            use_container_width=True,
                        )
                        st.download_button(
                            f"下载 {ontology} GO 气泡图 PNG",
                            data=bubbleplot_bytes,
                            file_name=f"GO_{ontology}_enrichment_bubbleplot.png",
                            mime="image/png",
                            key=f"download_go_{ontology}_bubbleplot",
                        )
                    st.caption(f"GO 气泡图生成用时：{time.perf_counter() - plot_started_at:.2f} 秒")

                st.download_button(
                    "下载全部结果 TSV",
                    data=results_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                    file_name="GO_enrichment_results.tsv",
                    mime="text/tab-separated-values",
                )
                st.download_button(
                    "下载显著结果 TSV",
                    data=sig_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                    file_name="GO_enrichment_results_sig.tsv",
                    mime="text/tab-separated-values",
                )
                st.download_button(
                    "下载分析摘要 TSV",
                    data=summary_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                    file_name="GO_enrichment_summary.tsv",
                    mime="text/tab-separated-values",
                )
            except Exception as e:
                st.error(f"GO 富集分析失败：{e}")
