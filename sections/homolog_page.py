import pandas as pd
import streamlit as st

from app_shared import EXAMPLE_CS_GENES, as_text, read_gene_ids, render_example_tools
from utils.db_query import (
    get_cs_self_all_hits,
    get_cs_self_best_hit,
    get_fielder_all_hits,
    get_fielder_best_hit,
)


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

    gene_ids = [g.strip() for g in read_gene_ids(uploaded_file, manual_input) if g.strip()]
    if not gene_ids:
        st.info("请上传文件或输入中国春基因号")
        st.stop()

    st.info(f"待查询基因数：{len(gene_ids)}")

    if st.button("开始查询同源基因", key="btn_blast"):
        progress = st.progress(0)
        status_text = st.empty()

        self_best_rows = []
        self_all_rows = []
        fld_best_rows = []
        fld_all_rows = []
        failed_genes = []

        for idx, input_gene_id in enumerate(gene_ids, 1):
            status_text.text(f"正在查询: {input_gene_id} ({idx}/{len(gene_ids)})")

            self_best_df = get_cs_self_best_hit(input_gene_id)
            self_all_df_one = get_cs_self_all_hits(input_gene_id)
            fld_best_df = get_fielder_best_hit(input_gene_id)
            fld_all_df_one = get_fielder_all_hits(input_gene_id)

            if self_best_df.empty and fld_best_df.empty:
                failed_genes.append(input_gene_id)

            if not self_best_df.empty:
                r = self_best_df.iloc[0]
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
                tmp = self_all_df_one.copy()
                tmp.insert(0, "input_gene_id", input_gene_id)
                self_all_rows.append(tmp)

            if not fld_best_df.empty:
                r = fld_best_df.iloc[0]
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
                tmp = fld_all_df_one.copy()
                tmp.insert(0, "input_gene_id", input_gene_id)
                fld_all_rows.append(tmp)

            progress.progress(idx / len(gene_ids))

        self_best_result_df = pd.DataFrame(
            self_best_rows,
            columns=["输入基因号", "中国春主键", "默认中国春自身同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"],
        )
        fld_best_result_df = pd.DataFrame(
            fld_best_rows,
            columns=["输入基因号", "中国春主键", "默认Fielder同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"],
        )

        self_all_df = pd.concat(self_all_rows, ignore_index=True) if self_all_rows else pd.DataFrame()
        fld_all_df = pd.concat(fld_all_rows, ignore_index=True) if fld_all_rows else pd.DataFrame()

        st.success("查询完成")
        st.subheader("中国春自身同源：默认最可信结果")
        st.dataframe(self_best_result_df, use_container_width=True)

        st.subheader("中国春自身同源：全部候选")
        if not self_all_df.empty:
            st.dataframe(self_all_df, use_container_width=True)
            st.download_button(
                "下载中国春自身同源全部候选 CSV",
                self_all_df.to_csv(index=False).encode("utf-8-sig"),
                "cs_self_homolog_all_hits.csv",
                "text/csv",
            )

        st.subheader("Fielder 同源：默认最可信结果")
        st.dataframe(fld_best_result_df, use_container_width=True)

        st.subheader("Fielder 同源：全部候选")
        if not fld_all_df.empty:
            st.dataframe(fld_all_df, use_container_width=True)
            st.download_button(
                "下载 Fielder 同源全部候选 CSV",
                fld_all_df.to_csv(index=False).encode("utf-8-sig"),
                "cs_to_fielder_all_hits.csv",
                "text/csv",
            )

        if failed_genes:
            st.warning(f"以下基因未找到同源结果: {', '.join(failed_genes)}")

        status_text.empty()
        progress.empty()
