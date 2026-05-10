import pandas as pd
import streamlit as st

from app_shared import EXAMPLE_CS_GENES, as_text, read_gene_ids, render_example_tools, show_large_input_notice
from utils.db_query import get_promoter


def render():
    st.header("中国春基因启动子抓取（本地数据库版）")

    uploaded_file = st.file_uploader("上传基因列表 TXT", type=["txt"], key="file_cs_promoter")
    render_example_tools(
        input_key="input_cs_promoter",
        example_text=as_text(EXAMPLE_CS_GENES),
        load_label="加载示例中国春基因",
        help_text="示例会自动填入下方输入框，用于测试中国春启动子提取。",
    )
    manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_cs_promoter")

    gene_ids = [g.strip() for g in read_gene_ids(uploaded_file, manual_input) if g.strip()]
    if not gene_ids:
        st.info("请上传文件或输入基因号")
        st.stop()

    st.info("当前数据库内置的是 promoter_2000，因此这里固定返回 ATG 上游 2000 bp 启动子序列。")
    st.info(f"待处理基因数: {len(gene_ids)}")
    show_large_input_notice(len(gene_ids), task_name="启动子序列提取", threshold=300)

    if st.button("开始抓取", key="btn_cs_promoter"):
        progress = st.progress(0)
        status_text = st.empty()

        fasta_records = []
        summary_rows = []
        failed_genes = []

        for idx, input_gene_id in enumerate(gene_ids, 1):
            status_text.text(f"正在查询: {input_gene_id} ({idx}/{len(gene_ids)})")

            promoter_df = get_promoter(input_gene_id)
            if promoter_df.empty:
                failed_genes.append(input_gene_id)
                summary_rows.append([input_gene_id, "", "", "", "", "", 0])
                progress.progress(idx / len(gene_ids))
                continue

            row = promoter_df.iloc[0]
            primary_id = row.get("primary_gene_id", "")
            chrom = row.get("chromosome", "")
            strand = row.get("strand", "")
            promoter_start = row.get("promoter_start", "")
            promoter_end = row.get("promoter_end", "")
            promoter_length = row.get("promoter_length", "")
            promoter_seq = row.get("promoter_sequence", "")

            fasta_records.append(
                f">{input_gene_id}|{primary_id}|promoter_2000|{chrom}:{promoter_start}-{promoter_end}|{strand}\n{promoter_seq}\n"
            )
            summary_rows.append([
                input_gene_id,
                primary_id if pd.notna(primary_id) else "",
                chrom if pd.notna(chrom) else "",
                promoter_start if pd.notna(promoter_start) else "",
                promoter_end if pd.notna(promoter_end) else "",
                strand if pd.notna(strand) else "",
                promoter_length if pd.notna(promoter_length) else 0,
            ])

            progress.progress(idx / len(gene_ids))

        summary_df = pd.DataFrame(
            summary_rows,
            columns=["输入基因号", "统一主键", "染色体", "启动子起点", "启动子终点", "链方向", "启动子长度"],
        )

        st.success("启动子抓取完成")
        st.dataframe(summary_df, use_container_width=True)
        if fasta_records:
            st.download_button("下载启动子 FASTA", "".join(fasta_records), "cs_promoter_sequences.fasta", "text/plain")
        else:
            st.warning("没有可下载的结果，请检查输入基因 ID 是否存在。")
        st.download_button(
            "下载启动子统计 CSV",
            summary_df.to_csv(index=False).encode("utf-8-sig"),
            "cs_promoter_summary.csv",
            "text/csv",
        )

        if failed_genes:
            st.warning(f"以下基因未获取到启动子序列: {', '.join(failed_genes)}")

        status_text.empty()
        progress.empty()
