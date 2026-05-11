import os

import streamlit as st

from app_shared import (
    get_go_example_text,
    get_go_mapping_paths,
    read_gene_ids,
    render_example_tools,
    run_cached_go_enrichment,
    show_large_input_notice,
)
from utils.go_enrichment import create_go_barplot_panel_bytes, create_go_bubbleplot_bytes


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

    gene_ids = [g.strip() for g in read_gene_ids(uploaded_file, manual_input) if g.strip()]
    if not gene_ids:
        st.info("请上传 DEG 文件或手动输入基因号")
        st.stop()
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
        st.caption("GO 富集图使用 qvalue 映射颜色；qvalue 越小颜色越偏红，越显著。气泡图横轴优先使用 RichFactor。")

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
                st.subheader("分析摘要")
                st.dataframe(summary_df, use_container_width=True)

                st.subheader("显著富集结果")
                if sig_df.empty:
                    st.warning("当前 FDR 阈值下没有显著富集条目，下面展示全部结果。")
                    st.dataframe(results_df, use_container_width=True)
                else:
                    st.dataframe(sig_df, use_container_width=True)

                plot_df = sig_df.copy()
                if plot_df.empty:
                    plot_df = results_df.copy()

                for ontology, label in [("BP", "BP"), ("CC", "CC"), ("MF", "MF")]:
                    if plot_df.empty or plot_df[plot_df["ontology"] == ontology].empty:
                        st.info(f"{label} 类别未检测到显著富集的 GO term。")

                show_bar = plot_type in {"Bar plot", "Both"}
                show_bubble = plot_type in {"Bubble plot", "Both"}

                if show_bar:
                    st.subheader("GO 富集条形图")
                    for ontology in ["BP", "CC", "MF"]:
                        barplot_bytes = create_go_barplot_panel_bytes(
                            df=plot_df,
                            ontology=ontology,
                            top_n=top_n,
                            qvalue_col="qvalue",
                            count_col="Count",
                            label_wrap_width=label_wrap_width,
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

                if show_bubble:
                    st.subheader("GO 富集气泡图")
                    for ontology in ["BP", "CC", "MF"]:
                        bubbleplot_bytes = create_go_bubbleplot_bytes(
                            df=plot_df,
                            ontology=ontology,
                            top_n=top_n,
                            qvalue_col="qvalue",
                            count_col="Count",
                            label_wrap_width=label_wrap_width,
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
