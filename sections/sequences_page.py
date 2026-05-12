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
from utils.db_query import get_gene_sequence_resources_many, resolve_primary_gene_ids_many


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

    gene_ids, cleanup_info = preprocess_gene_ids(uploaded_file, manual_input)
    if not gene_ids:
        st.info("请上传文件或输入基因号")
        st.stop()
    show_input_cleanup_notice(cleanup_info)

    st.info(f"待查询基因数：{len(gene_ids)}")
    show_large_input_notice(len(gene_ids), task_name="序列查询", threshold=300)

    if st.button("获取序列（cDNA / CDS / Protein）", key="btn_sequences"):
        started_at = time.perf_counter()
        cdna_records = []
        cds_records = []
        protein_records = []
        failed_genes = []
        summary_rows = []

        with st.spinner("正在批量获取 cDNA / CDS / protein 序列，请稍候..."):
            id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
            all_seq_df = get_gene_sequence_resources_many(tuple(gene_ids))

        grouped = {
            input_gene_id: group.copy()
            for input_gene_id, group in all_seq_df.groupby("input_gene_id", sort=False)
        } if not all_seq_df.empty else {}

        for input_gene_id in gene_ids:
            primary_gene_id = id_map.get(input_gene_id)
            seq_df = grouped.get(input_gene_id, pd.DataFrame())
            if seq_df.empty:
                failed_genes.append(input_gene_id)
                summary_rows.append([input_gene_id, "", 0, 0, 0])
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

        summary_df = pd.DataFrame(
            summary_rows,
            columns=["输入基因号", "统一主键", "cDNA条数", "CDS条数", "Protein条数"],
        )

        st.success("序列查询完成")
        st.caption(f"查询用时：{time.perf_counter() - started_at:.2f} 秒")
        show_dataframe_preview(summary_df, label="序列统计结果", key="show_all_sequence_summary")
        if cdna_records:
            st.download_button("下载 cDNA FASTA", "".join(cdna_records), "cdna_sequences.fasta", "text/plain")
        if cds_records:
            st.download_button("下载 CDS FASTA", "".join(cds_records), "cds_sequences.fasta", "text/plain")
        if protein_records:
            st.download_button("下载 Protein FASTA", "".join(protein_records), "protein_sequences.fasta", "text/plain")
        if not (cdna_records or cds_records or protein_records):
            st.warning("没有可下载的结果，请检查输入基因 ID 是否存在。")
        st.download_button(
            "下载序列统计 CSV",
            summary_df.to_csv(index=False).encode("utf-8-sig"),
            "sequence_summary.csv",
            "text/csv",
        )

        if failed_genes:
            st.warning(f"以下基因未获取到序列: {', '.join(failed_genes)}")
