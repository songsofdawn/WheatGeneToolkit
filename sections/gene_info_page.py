import time

import pandas as pd
import streamlit as st

from app_shared import (
    EXAMPLE_CS_GENES,
    as_text,
    preprocess_gene_ids,
    render_example_tools,
    show_dataframe_preview,
    show_input_cleanup_notice,
    show_large_input_notice,
)
from utils.db_query import get_gene_annotations_many, get_gene_core_many, resolve_primary_gene_ids_many


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

    gene_ids, cleanup_info = preprocess_gene_ids(uploaded_file, manual_input)
    if not gene_ids:
        st.info("请上传文件或输入基因号")
        st.stop()
    show_input_cleanup_notice(cleanup_info)

    st.info(f"待查询基因数：{len(gene_ids)}")
    show_large_input_notice(len(gene_ids), task_name="基因功能注释查询", threshold=500)

    if st.button("开始查询", key="btn_gene_info"):
        started_at = time.perf_counter()
        results = []

        with st.spinner("正在批量查询基因功能注释，请稍候..."):
            id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
            anno_df = get_gene_annotations_many(tuple(gene_ids))
            core_df = get_gene_core_many(tuple(gene_ids))

        anno_by_input = {
            row["input_gene_id"]: row
            for _, row in anno_df.drop_duplicates("input_gene_id").iterrows()
        } if not anno_df.empty else {}
        core_by_input = {
            row["input_gene_id"]: row
            for _, row in core_df.drop_duplicates("input_gene_id").iterrows()
        } if not core_df.empty else {}

        missing_genes = []
        for input_gene_id in gene_ids:
            primary_gene_id = id_map.get(input_gene_id)
            anno_row = anno_by_input.get(input_gene_id)
            core_row = core_by_input.get(input_gene_id)

            third_id = core_row.get("gene_id_v3", None) if core_row is not None else None
            description_en = anno_row.get("description_en", None) if anno_row is not None else None
            description_zh = anno_row.get("description_zh", None) if anno_row is not None else None

            if primary_gene_id is None:
                missing_genes.append(input_gene_id)
            results.append([
                input_gene_id,
                primary_gene_id if primary_gene_id is not None else "未找到",
                third_id if pd.notna(third_id) else "",
                description_en if pd.notna(description_en) else "",
                description_zh if pd.notna(description_zh) else "",
            ])

        df = pd.DataFrame(
            results,
            columns=["输入基因号", "统一主键", "三代基因号", "功能描述（英文）", "功能描述（中文）"],
        )

        st.success("查询完成")
        st.caption(f"查询用时：{time.perf_counter() - started_at:.2f} 秒")
        show_dataframe_preview(df, label="基因功能注释结果", key="show_all_gene_info")
        st.download_button(
            "下载结果 CSV",
            df.to_csv(index=False).encode("utf-8-sig"),
            "gene_info_from_sqlite.csv",
            "text/csv",
        )
        if missing_genes:
            st.warning(f"以下基因未找到: {', '.join(missing_genes)}")
