import os

import streamlit as st

from app_shared import (
    get_kegg_example_text,
    get_kegg_mapping_paths,
    read_gene_ids,
    render_example_tools,
    run_cached_kegg_enrichment,
    show_large_input_notice,
)
from utils.kegg_enrichment import create_kegg_barplot_bytes, create_kegg_bubbleplot_bytes


def render():
    st.header("KEGG 富集分析")
    st.caption("输入 DEG 列表（一行一个基因号），基于本地 gene-KO 和 KEGG KO-pathway 映射进行 KEGG 富集分析。")

    uploaded_file = st.file_uploader(
        "上传 DEG TXT 文件（一行一个基因号）",
        type=["txt"],
        key="file_kegg_deg",
    )
    render_example_tools(
        input_key="input_kegg_deg",
        example_text=get_kegg_example_text(),
        load_label="加载 KEGG 示例 DEG",
        download_label="下载 KEGG 示例 DEG 文件",
        file_name="example_kegg_deg_genes.txt",
        help_text="示例 DEG 从 data/example_kegg_deg_genes.txt 读取；你可以直接修改该文件来更新 KEGG 示例。",
    )
    manual_input = st.text_area("或者手动输入 DEG（每行一个）", height=200, key="input_kegg_deg")

    gene_ids = [g.strip() for g in read_gene_ids(uploaded_file, manual_input) if g.strip()]
    if not gene_ids:
        st.info("请上传 DEG 文件或手动输入基因号")
        st.stop()

    st.info(f"待分析基因数: {len(gene_ids)}")
    show_large_input_notice(len(gene_ids), task_name="KEGG 富集分析", threshold=2000)

    col1, col2, col3 = st.columns(3)
    with col1:
        top_n = st.number_input("图中显示前 N 个通路", min_value=3, max_value=50, value=20, step=1)
    with col2:
        pvalue_cutoff = st.number_input("P-value 阈值", min_value=0.0001, max_value=1.0, value=0.05, step=0.01, format="%.4f")
    with col3:
        use_sig_only = st.checkbox(
            "优先显示显著通路",
            value=True,
            help="如果存在 p-value 小于阈值的通路，则优先使用显著通路绘图；否则展示全部结果中 p-value 最小的通路。",
        )

    col4, col5 = st.columns(2)
    with col4:
        min_size = st.number_input("最小 pathway KO 数", min_value=1, max_value=100, value=3, step=1)
    with col5:
        max_size = st.number_input("最大 pathway KO 数", min_value=10, max_value=5000, value=500, step=10)

    with st.expander("高级绘图参数"):
        col6, col7, col8 = st.columns(3)
        with col6:
            clip_minus_log10_p = st.checkbox(
                "截断极端 -log₁₀(P-value)",
                value=True,
                help="建议开启，避免极端小 p-value 把横轴拉得过长。",
            )
        with col7:
            clip_quantile = st.number_input("截断分位数", min_value=0.50, max_value=1.00, value=0.95, step=0.01, format="%.2f")
        with col8:
            label_wrap_width = st.number_input("通路名称换行宽度", min_value=20, max_value=80, value=42, step=2)

        col9, col10 = st.columns(2)
        with col9:
            bubble_min_size = st.number_input("最小气泡大小", min_value=5, max_value=200, value=25, step=5)
        with col10:
            bubble_max_size = st.number_input("最大气泡大小", min_value=50, max_value=800, value=220, step=10)

    kegg_paths = get_kegg_mapping_paths()
    gene2ko_path = kegg_paths["gene2ko"]
    ko2pathway_path = kegg_paths["ko2pathway"]
    pathway2name_path = kegg_paths["pathway2name"]

    required_files = [gene2ko_path, ko2pathway_path, pathway2name_path]
    missing_files = [fp for fp in required_files if not os.path.exists(fp)]
    if missing_files:
        st.error("以下 KEGG 注释文件不存在，请检查 data/kegg_mapping/ 目录：")
        for fp in missing_files:
            st.code(fp)
        st.stop()

    with st.expander("当前使用的 KEGG 本地注释文件"):
        st.code(gene2ko_path)
        st.code(ko2pathway_path)
        st.code(pathway2name_path)

    if st.button("开始 KEGG 富集分析", key="btn_kegg_enrichment"):
        with st.spinner("正在进行 KEGG 富集分析，请稍候..."):
            try:
                results_df, sig_df, summary_df = run_cached_kegg_enrichment(
                    gene_ids=tuple(gene_ids),
                    gene2ko_path=gene2ko_path,
                    ko2pathway_path=ko2pathway_path,
                    pathway2name_path=pathway2name_path,
                    min_size=min_size,
                    max_size=max_size,
                    pvalue_cutoff=pvalue_cutoff,
                )

                if results_df.empty:
                    st.warning("没有得到任何 KEGG 富集结果，请检查输入基因 ID 是否与 gene2ko_clean.tsv 中的基因 ID 一致。")
                    st.subheader("分析摘要")
                    st.dataframe(summary_df, use_container_width=True)
                    st.stop()

                st.success("KEGG 富集分析完成")
                st.subheader("分析摘要")
                st.dataframe(summary_df, use_container_width=True)

                st.subheader("显著富集结果")
                if sig_df.empty:
                    st.warning("当前 P-value 阈值下没有显著富集通路，下面展示全部结果。")
                    st.dataframe(results_df, use_container_width=True)
                else:
                    st.dataframe(sig_df, use_container_width=True)

                plot_df = sig_df if use_sig_only and not sig_df.empty else results_df
                bubble_bytes = create_kegg_bubbleplot_bytes(
                    df=plot_df,
                    top_n=top_n,
                    bubble_min_size=bubble_min_size,
                    bubble_max_size=bubble_max_size,
                    clip_minus_log10_p=clip_minus_log10_p,
                    clip_mode="quantile",
                    clip_quantile=clip_quantile,
                    clip_fixed_value=30,
                    min_x_cap=5,
                    label_wrap_width=label_wrap_width,
                )
                barplot_bytes = create_kegg_barplot_bytes(
                    df=plot_df,
                    top_n=top_n,
                    clip_minus_log10_p=clip_minus_log10_p,
                    clip_mode="quantile",
                    clip_quantile=clip_quantile,
                    clip_fixed_value=30,
                    min_x_cap=5,
                    label_wrap_width=label_wrap_width,
                )

                if bubble_bytes is not None:
                    st.subheader("KEGG 富集气泡图")
                    st.image(bubble_bytes, caption="KEGG enrichment bubble plot", use_container_width=True)
                    st.download_button(
                        "下载 KEGG 气泡图 PNG",
                        data=bubble_bytes,
                        file_name="KEGG_enrichment_bubbleplot.png",
                        mime="image/png",
                    )

                if barplot_bytes is not None:
                    st.subheader("KEGG 富集条形图")
                    st.image(barplot_bytes, caption="KEGG enrichment barplot", use_container_width=True)
                    st.download_button(
                        "下载 KEGG 条形图 PNG",
                        data=barplot_bytes,
                        file_name="KEGG_enrichment_barplot.png",
                        mime="image/png",
                    )

                st.download_button(
                    "下载全部 KEGG 富集结果 TSV",
                    data=results_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                    file_name="KEGG_enrichment_results_all.tsv",
                    mime="text/tab-separated-values",
                )
                st.download_button(
                    "下载显著 KEGG 富集结果 TSV",
                    data=sig_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                    file_name="KEGG_enrichment_results_sig.tsv",
                    mime="text/tab-separated-values",
                )
                st.download_button(
                    "下载 KEGG 分析摘要 TSV",
                    data=summary_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                    file_name="KEGG_enrichment_summary.tsv",
                    mime="text/tab-separated-values",
                )
            except Exception as e:
                st.error(f"KEGG 富集分析失败：{e}")
