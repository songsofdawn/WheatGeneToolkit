import pandas as pd
import streamlit as st

from app_shared import EXAMPLE_CS_GENES, as_text, read_gene_ids, render_example_tools
from utils.db_query import get_primary_gene_id, get_sequences


def render():
    st.header("cDNA / CDS / Protein 下载（本地数据库版）")

    uploaded_file = st.file_uploader("上传 TXT 文件（基因号）", type=["txt"], key="file_sequences")
    render_example_tools(
        input_key="input_sequences",
        example_text=as_text(EXAMPLE_CS_GENES),
        load_label="加载示例中国春基因",
        help_text="示例会自动填入下方输入框，可直接获取 cDNA / CDS / protein。",
    )
    manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_sequences")

    gene_ids = [g.strip() for g in read_gene_ids(uploaded_file, manual_input) if g.strip()]
    if not gene_ids:
        st.info("请上传文件或输入基因号")
        st.stop()

    st.info(f"待查询基因数：{len(gene_ids)}")

    if st.button("获取序列（cDNA / CDS / Protein）", key="btn_sequences"):
        cdna_records = []
        cds_records = []
        protein_records = []
        failed_genes = []
        summary_rows = []

        progress = st.progress(0)
        status_text = st.empty()

        for idx, input_gene_id in enumerate(gene_ids, 1):
            status_text.text(f"正在处理: {input_gene_id} ({idx}/{len(gene_ids)})")

            primary_gene_id = get_primary_gene_id(input_gene_id)
            seq_df = get_sequences(input_gene_id)

            if seq_df.empty:
                failed_genes.append(input_gene_id)
                summary_rows.append([input_gene_id, "", 0, 0, 0])
                progress.progress(idx / len(gene_ids))
                continue

            cdna_df = seq_df[seq_df["sequence_type"] == "cdna"].copy()
            cds_df = seq_df[seq_df["sequence_type"] == "cds"].copy()
            protein_df = seq_df[seq_df["sequence_type"] == "protein"].copy()

            summary_rows.append([
                input_gene_id,
                primary_gene_id if primary_gene_id else "",
                len(cdna_df),
                len(cds_df),
                len(protein_df),
            ])

            for _, row in cdna_df.iterrows():
                tx_id = row.get("transcript_id", "")
                seq = row.get("sequence", "")
                cdna_records.append(f">{input_gene_id}|{primary_gene_id}|{tx_id}|cdna\n{seq}\n")

            for _, row in cds_df.iterrows():
                tx_id = row.get("transcript_id", "")
                seq = row.get("sequence", "")
                cds_records.append(f">{input_gene_id}|{primary_gene_id}|{tx_id}|cds\n{seq}\n")

            for _, row in protein_df.iterrows():
                tx_id = row.get("transcript_id", "")
                seq = row.get("sequence", "")
                protein_records.append(f">{input_gene_id}|{primary_gene_id}|{tx_id}|protein\n{seq}\n")

            if len(cdna_df) == 0 and len(cds_df) == 0 and len(protein_df) == 0:
                failed_genes.append(input_gene_id)

            progress.progress(idx / len(gene_ids))

        summary_df = pd.DataFrame(
            summary_rows,
            columns=["输入基因号", "统一主键", "cDNA条数", "CDS条数", "Protein条数"],
        )

        st.success("序列查询完成")
        st.dataframe(summary_df, use_container_width=True)
        st.download_button("下载 cDNA FASTA", "".join(cdna_records), "cdna_sequences.fasta", "text/plain")
        st.download_button("下载 CDS FASTA", "".join(cds_records), "cds_sequences.fasta", "text/plain")
        st.download_button("下载 Protein FASTA", "".join(protein_records), "protein_sequences.fasta", "text/plain")
        st.download_button(
            "下载序列统计 CSV",
            summary_df.to_csv(index=False).encode("utf-8-sig"),
            "sequence_summary.csv",
            "text/csv",
        )

        if failed_genes:
            st.warning(f"以下基因未获取到序列: {', '.join(failed_genes)}")

        status_text.empty()
        progress.empty()
