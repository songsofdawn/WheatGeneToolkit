import pandas as pd
import streamlit as st

from app_shared import EXAMPLE_CS_GENES, as_text, read_gene_ids, render_example_tools, show_large_input_notice
from utils.db_query import get_gene_annotation, get_gene_core, get_primary_gene_id


def render():
    st.header("基因功能注释及三代基因号转换（本地数据库版）")

    uploaded_file = st.file_uploader(
        "上传 TXT 文件（基因号一行一个，不要加 .1）",
        type=["txt"],
        key="file_gene_info",
    )
    render_example_tools(
        input_key="input_gene_info",
        example_text=as_text(EXAMPLE_CS_GENES),
        load_label="加载示例中国春基因",
        help_text="示例会自动填入下方输入框，可直接点击开始查询。",
    )
    manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_gene_info")

    gene_ids = [g.strip() for g in read_gene_ids(uploaded_file, manual_input) if g.strip()]
    if not gene_ids:
        st.info("请上传文件或输入基因号")
        st.stop()

    st.info(f"待查询基因数：{len(gene_ids)}")
    show_large_input_notice(len(gene_ids), task_name="基因功能注释查询", threshold=500)

    if st.button("开始查询", key="btn_gene_info"):
        progress_bar = st.progress(0)
        status_text = st.empty()
        results = []

        for idx, input_gene_id in enumerate(gene_ids, 1):
            status_text.text(f"正在查询：{input_gene_id} ({idx}/{len(gene_ids)})")

            primary_gene_id = get_primary_gene_id(input_gene_id)
            anno_df = get_gene_annotation(input_gene_id)
            core_df = get_gene_core(input_gene_id)

            third_id = None
            description_en = None
            description_zh = None

            if not anno_df.empty:
                description_en = anno_df.iloc[0].get("description_en", None)
                description_zh = anno_df.iloc[0].get("description_zh", None)

            if not core_df.empty:
                third_id = core_df.iloc[0].get("gene_id_v3", None)

            results.append([
                input_gene_id,
                primary_gene_id if primary_gene_id is not None else "未找到",
                third_id if pd.notna(third_id) else "",
                description_en if pd.notna(description_en) else "",
                description_zh if pd.notna(description_zh) else "",
            ])

            progress_bar.progress(idx / len(gene_ids))

        df = pd.DataFrame(
            results,
            columns=["输入基因号", "统一主键", "三代基因号", "功能描述（英文）", "功能描述（中文）"],
        )

        st.success("查询完成")
        st.dataframe(df, use_container_width=True)
        st.download_button(
            "下载结果 CSV",
            df.to_csv(index=False).encode("utf-8-sig"),
            "gene_info_from_sqlite.csv",
            "text/csv",
        )

        status_text.empty()
        progress_bar.empty()
