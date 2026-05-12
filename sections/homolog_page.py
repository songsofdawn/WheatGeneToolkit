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
from utils.db_query import get_homologs_many


def render():
    st.header("中国春同源基因检索（自身同源 + Fielder）")

    uploaded_file = st.file_uploader("上传 TXT 文件（基因号一行一个）", type=["txt"], key="file_blast")
    render_example_tools(
        input_key="input_blast",
        example_text=as_text(EXAMPLE_CS_GENES),
        load_label="加载示例中国春基因",
        help_text="示例会自动填入下方输入框，可直接测试自身同源和 Fielder 同源查询。",
    )
    manual_input = st.text_area("或者手动输入中国春基因号（每行一个）", height=200, key="input_blast")

    gene_ids, cleanup_info = preprocess_gene_ids(uploaded_file, manual_input)
    if not gene_ids:
        st.info("请上传文件或输入中国春基因号")
        st.stop()
    show_input_cleanup_notice(cleanup_info)

    st.info(f"待查询基因数：{len(gene_ids)}")
    show_large_input_notice(len(gene_ids), task_name="同源基因查询", threshold=500)

    if st.button("开始查询同源基因", key="btn_blast"):
        started_at = time.perf_counter()

        self_best_rows = []
        fld_best_rows = []

        with st.spinner("正在批量查询同源基因，请稍候..."):
            homolog_results = get_homologs_many(tuple(gene_ids))

        self_best_df = homolog_results["self_best"]
        fld_best_df = homolog_results["fielder_best"]
        self_all_df = homolog_results["self_all"]
        fld_all_df = homolog_results["fielder_all"]
        failed_genes = homolog_results["missing"]

        self_best_by_input = {
            input_gene_id: group.iloc[0]
            for input_gene_id, group in self_best_df.groupby("input_gene_id", sort=False)
        } if not self_best_df.empty else {}
        fld_best_by_input = {
            input_gene_id: group.iloc[0]
            for input_gene_id, group in fld_best_df.groupby("input_gene_id", sort=False)
        } if not fld_best_df.empty else {}

        for input_gene_id in gene_ids:
            if input_gene_id in self_best_by_input:
                r = self_best_by_input[input_gene_id]
                self_best_rows.append([
                    input_gene_id,
                    r.get("cs_gene_id", ""),
                    r.get("self_homolog_gene_id", ""),
                    r.get("homolog_type", ""),
                    r.get("confidence", ""),
                    r.get("same_subgenome", ""),
                    r.get("same_chr_num", ""),
                    r.get("priority_score", ""),
                ])

            if input_gene_id in fld_best_by_input:
                r = fld_best_by_input[input_gene_id]
                fld_best_rows.append([
                    input_gene_id,
                    r.get("cs_gene_id", ""),
                    r.get("fielder_gene_id", ""),
                    r.get("homolog_type", ""),
                    r.get("confidence", ""),
                    r.get("same_subgenome", ""),
                    r.get("same_chr_num", ""),
                    r.get("priority_score", ""),
                ])

        self_best_result_df = pd.DataFrame(
            self_best_rows,
            columns=["输入基因号", "中国春主键", "默认中国春自身同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"],
        )
        fld_best_result_df = pd.DataFrame(
            fld_best_rows,
            columns=["输入基因号", "中国春主键", "默认Fielder同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"],
        )

        st.success("查询完成")
        st.caption(f"查询用时：{time.perf_counter() - started_at:.2f} 秒")
        st.subheader("中国春自身同源：默认最可信结果")
        show_dataframe_preview(self_best_result_df, label="中国春自身同源默认结果", key="show_all_self_best")

        st.subheader("中国春自身同源：全部候选")
        if not self_all_df.empty:
            show_dataframe_preview(self_all_df, label="中国春自身同源全部候选", key="show_all_self_all")
            st.download_button(
                "下载中国春自身同源全部候选 CSV",
                self_all_df.to_csv(index=False).encode("utf-8-sig"),
                "cs_self_homolog_all_hits.csv",
                "text/csv",
            )

        st.subheader("Fielder 同源：默认最可信结果")
        show_dataframe_preview(fld_best_result_df, label="Fielder 同源默认结果", key="show_all_fielder_best")

        st.subheader("Fielder 同源：全部候选")
        if not fld_all_df.empty:
            show_dataframe_preview(fld_all_df, label="Fielder 同源全部候选", key="show_all_fielder_all")
            st.download_button(
                "下载 Fielder 同源全部候选 CSV",
                fld_all_df.to_csv(index=False).encode("utf-8-sig"),
                "cs_to_fielder_all_hits.csv",
                "text/csv",
            )

        if failed_genes:
            st.warning(f"以下基因未找到同源结果: {', '.join(failed_genes)}")
