import os

import pandas as pd
import streamlit as st

from app_shared import DATA_DIR
from utils.volcano_plot import (
    auto_detect_columns,
    figure_to_png_bytes,
    infer_auto_y_cap,
    plot_volcano,
    prepare_volcano_data,
    read_volcano_table,
)


EXAMPLE_VOLCANO_PATH = os.path.join(DATA_DIR, "example_volcano_plot_genes_list.txt")
EXAMPLE_VOLCANO_REL = "data/example_volcano_plot_genes_list.txt"


def _load_example_text() -> str:
    if not os.path.exists(EXAMPLE_VOLCANO_PATH):
        return ""
    with open(EXAMPLE_VOLCANO_PATH, "r", encoding="utf-8") as f:
        return f.read()


def _select_index(columns, detected_col):
    if detected_col in columns:
        return columns.index(detected_col)
    return 0


def _optional_select_index(options, detected_col):
    if detected_col in options:
        return options.index(detected_col)
    return 0


def _txt_download(series: pd.Series) -> bytes:
    values = [str(value) for value in series.dropna().astype(str).tolist() if str(value).strip()]
    return ("\n".join(values) + ("\n" if values else "")).encode("utf-8")


def render():
    st.header("火山图分析 Volcano Plot")
    st.caption("本模块用于根据差异表达结果表绘制火山图。输入表格只需要包含 log₂FC 和 p-value 两列；gene_id 列可选，用于下载基因列表和标注。")

    uploaded_file = st.file_uploader(
        "上传差异表达结果文件（txt / tsv / csv）",
        type=["txt", "tsv", "csv"],
        key="file_volcano",
    )

    if st.button("加载示例：MYC2 48H 差异表达火山图数据", key="load_volcano_example"):
        example_text = _load_example_text()
        if example_text:
            st.session_state["input_volcano_table"] = example_text
            st.session_state["volcano_example_loaded"] = True
        else:
            st.error(f"未找到示例文件：{EXAMPLE_VOLCANO_REL}")

    if st.session_state.get("volcano_example_loaded"):
        st.success(f"已加载示例文件：{EXAMPLE_VOLCANO_REL}")

    manual_input = st.text_area(
        "或者粘贴差异表达结果表（第一行为表头）",
        height=220,
        key="input_volcano_table",
    )

    if uploaded_file is None and not manual_input.strip():
        st.info("请上传文件、粘贴表格文本，或点击加载示例数据。")
        st.stop()

    try:
        source = uploaded_file if uploaded_file is not None else manual_input
        df = read_volcano_table(source)
    except Exception as exc:
        st.error(f"读取表格失败：{exc}")
        st.stop()

    if df.empty:
        st.error("读取到的表格为空，请检查输入内容。")
        st.stop()

    st.subheader("数据预览")
    st.dataframe(df.head(20), use_container_width=True)
    st.caption(f"原始行数：{len(df)}；原始列数：{len(df.columns)}")

    detected = auto_detect_columns(df)
    columns = list(df.columns)
    optional_gene_options = ["不使用 gene_id 列"] + columns

    st.subheader("列选择")
    col1, col2, col3 = st.columns(3)
    with col1:
        gene_id_col = st.selectbox(
            "gene_id 列（可选）",
            optional_gene_options,
            index=_optional_select_index(optional_gene_options, detected["gene_id_col"]),
        )
        if gene_id_col == "不使用 gene_id 列":
            gene_id_col = None
    with col2:
        log2fc_col = st.selectbox(
            "log₂FC 列",
            columns,
            index=_select_index(columns, detected["log2fc_col"]),
        )
    with col3:
        pvalue_col = st.selectbox(
            "p-value 列",
            columns,
            index=_select_index(columns, detected["pvalue_col"]),
        )

    if not all([log2fc_col, pvalue_col]):
        st.warning("列识别不完整，请手动选择 log₂FC 和 p-value 列。gene_id 列可以不选。")

    st.subheader("绘图参数")
    col4, col5, col6 = st.columns(3)
    with col4:
        log2fc_cutoff = st.number_input("log₂FC cutoff", min_value=0.1, max_value=5.0, value=1.0, step=0.1)
    with col5:
        pvalue_cutoff = st.number_input("p-value cutoff", min_value=1e-10, max_value=1.0, value=0.05, step=0.01, format="%.6f")
    with col6:
        top_label_n = st.number_input("top label N", min_value=0, max_value=50, value=0, step=1)

    col7, col8, col9 = st.columns(3)
    with col7:
        point_size = st.number_input("point size", min_value=4, max_value=80, value=14, step=1)
    with col8:
        alpha = st.number_input("alpha", min_value=0.1, max_value=1.0, value=0.75, step=0.05, format="%.2f")
    with col9:
        cap_mode = st.selectbox("y-axis cap", ["自动", "20", "50", "100", "不截断"], index=0)

    st.caption("如果存在极小 p-value，y 轴可能被少数极端点拉高。默认会使用软压缩坐标轴保留主体点云，同时让极端点仍显示在图内。")

    title = st.text_input("图标题", value="Volcano Plot")

    if st.button("生成火山图", key="btn_volcano_plot"):
        if cap_mode == "自动":
            max_neg_log10_pvalue = None
        elif cap_mode == "不截断":
            max_neg_log10_pvalue = float("inf")
        else:
            max_neg_log10_pvalue = float(cap_mode)

        volcano_df_no_cap = prepare_volcano_data(
            df,
            gene_id_col=gene_id_col,
            log2fc_col=log2fc_col,
            pvalue_col=pvalue_col,
            log2fc_cutoff=log2fc_cutoff,
            pvalue_cutoff=pvalue_cutoff,
            max_neg_log10_pvalue=float("inf"),
        )

        if cap_mode == "自动" and not volcano_df_no_cap.empty:
            max_neg_log10_pvalue = infer_auto_y_cap(volcano_df_no_cap["neg_log10_pvalue_raw"])

        volcano_df = prepare_volcano_data(
            df,
            gene_id_col=gene_id_col,
            log2fc_col=log2fc_col,
            pvalue_col=pvalue_col,
            log2fc_cutoff=log2fc_cutoff,
            pvalue_cutoff=pvalue_cutoff,
            max_neg_log10_pvalue=max_neg_log10_pvalue,
        )

        if volcano_df.empty:
            st.error("清洗后没有有效数据。请检查 log₂FC / p-value 列是否为数字，且 p-value 是否在 (0, 1] 范围内。")
            st.stop()

        filtered_count = len(df) - len(volcano_df)
        up_count = int((volcano_df["regulation"] == "Up-regulated").sum())
        down_count = int((volcano_df["regulation"] == "Down-regulated").sum())
        not_sig_count = int((volcano_df["regulation"] == "Not significant").sum())
        capped_count = int(volcano_df["is_y_capped"].sum())

        st.subheader("统计摘要")
        m1, m2, m3, m4, m5, m6 = st.columns(6)
        m1.metric("原始行数", len(df))
        m2.metric("有效基因数", len(volcano_df))
        m3.metric("Up-regulated", up_count)
        m4.metric("Down-regulated", down_count)
        m5.metric("Not significant", not_sig_count)
        m6.metric("过滤行数", filtered_count)

        if capped_count:
            st.info(f"有 {capped_count} 个点被 y-axis cap 截断显示；这些点的真实 -log₁₀(p-value) 更高。")

        try:
            fig = plot_volcano(
                volcano_df,
                log2fc_cutoff=log2fc_cutoff,
                pvalue_cutoff=pvalue_cutoff,
                title=title,
                top_label_n=top_label_n,
                max_neg_log10_pvalue=max_neg_log10_pvalue,
                point_size=point_size,
                alpha=alpha,
            )
        except Exception as exc:
            st.error(f"绘图失败：{exc}")
            st.stop()

        png_bytes = figure_to_png_bytes(fig)
        st.subheader("火山图")
        st.image(png_bytes, caption="Volcano Plot", use_container_width=True)
        st.download_button(
            "下载火山图 PNG",
            data=png_bytes,
            file_name="volcano_plot.png",
            mime="image/png",
        )

        st.subheader("清洗后的结果表")
        st.dataframe(volcano_df, use_container_width=True)

        up_genes = volcano_df.loc[volcano_df["regulation"] == "Up-regulated", "gene_id"]
        down_genes = volcano_df.loc[volcano_df["regulation"] == "Down-regulated", "gene_id"]
        sig_genes = volcano_df.loc[volcano_df["regulation"].isin(["Up-regulated", "Down-regulated"]), "gene_id"]

        st.download_button(
            "下载清洗后的 volcano result CSV",
            data=volcano_df.to_csv(index=False).encode("utf-8-sig"),
            file_name="volcano_result_cleaned.csv",
            mime="text/csv",
        )
        st.download_button(
            "下载上调基因列表 txt",
            data=_txt_download(up_genes),
            file_name="volcano_upregulated_genes.txt",
            mime="text/plain",
        )
        st.download_button(
            "下载下调基因列表 txt",
            data=_txt_download(down_genes),
            file_name="volcano_downregulated_genes.txt",
            mime="text/plain",
        )
        st.download_button(
            "下载显著差异基因列表 txt",
            data=_txt_download(sig_genes),
            file_name="volcano_significant_genes.txt",
            mime="text/plain",
        )
