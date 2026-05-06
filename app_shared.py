import os

import streamlit as st

from utils.go_enrichment import run_go_enrichment
from utils.kegg_enrichment import run_kegg_enrichment

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
GO_MAPPING_DIR = os.path.join(DATA_DIR, "go_mapping")
KEGG_MAPPING_DIR = os.path.join(DATA_DIR, "kegg_mapping")
TIP_CODE_PATH = os.path.join(DATA_DIR, "TipCode.jpg")
README_PATH = os.path.join(BASE_DIR, "README.md")
EXAMPLE_GO_DEG_PATH = os.path.join(DATA_DIR, "example_go_deg_genes.txt")
EXAMPLE_KEGG_DEG_PATH = os.path.join(DATA_DIR, "example_kegg_deg_genes.txt")

EXAMPLE_CS_GENES = [
    "TraesCS2D02G571200",
    "TraesCS5B02G233300",
    "TraesCS1A02G000100",
    "TraesCS3D02G114000",
    "TraesCS4A02G109200",
]

EXAMPLE_FIELDER_GENES = [
    "TraesFLD5B01G105200",
    "TraesFLD1A01G000100",
    "TraesFLD2B01G000500",
]

DEFAULT_GO_DEG_GENES = [
    "TraesCS2D02G571200",
    "TraesCS5B02G233300",
    "TraesCS1A02G000100",
    "TraesCS3D02G114000",
    "TraesCS4A02G109200",
    "TraesCS2A02G286700",
    "TraesCS4B02G195000",
    "TraesCS4D02G195700",
    "TraesCS5A02G212800",
    "TraesCS5B02G211000",
    "TraesCS7A02G201100",
    "TraesCS7A02G201200",
    "TraesCS7A02G201300",
    "TraesCS7B02G107700",
    "TraesCS7B02G107800",
]

DEFAULT_KEGG_DEG_GENES = [
    "TraesCS5B02G233300",
    "TraesCS3D02G114000",
    "TraesCS2D02G571200",
    "TraesCS1A02G000100",
    "TraesCS4A02G109200",
    "TraesCS2A02G286700",
    "TraesCS4A02G007800",
    "TraesCS4B02G195000",
    "TraesCS4D02G195700",
    "TraesCS5A02G212800",
    "TraesCS5A02G533000",
    "TraesCS5A02G533100",
    "TraesCS5B02G211000",
    "TraesCS7A02G201100",
    "TraesCS7A02G201200",
    "TraesCS7A02G201300",
    "TraesCS7A02G201400",
    "TraesCS7A02G201600",
]


def read_gene_ids(uploaded_file, manual_input):
    if uploaded_file:
        return uploaded_file.read().decode("utf-8").splitlines()
    if manual_input.strip():
        return manual_input.strip().splitlines()
    return []


def as_text(items):
    return "\n".join(items)


def load_example_gene_text(file_path, fallback_items=None):
    if os.path.exists(file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            genes = []
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                genes.append(line)
        return as_text(genes)

    if fallback_items is not None:
        return as_text(fallback_items)

    return ""


def get_go_example_text():
    return load_example_gene_text(EXAMPLE_GO_DEG_PATH, DEFAULT_GO_DEG_GENES)


def get_kegg_example_text():
    return load_example_gene_text(EXAMPLE_KEGG_DEG_PATH, DEFAULT_KEGG_DEG_GENES)


def render_example_tools(
    input_key,
    example_text,
    load_label="加载示例数据",
    download_label=None,
    file_name=None,
    help_text=None,
):
    if download_label and file_name:
        col1, col2 = st.columns([1, 1])
        with col1:
            if st.button(load_label, key=f"load_example_{input_key}"):
                st.session_state[input_key] = example_text
        with col2:
            st.download_button(
                download_label,
                data=example_text.encode("utf-8"),
                file_name=file_name,
                mime="text/plain",
                key=f"download_example_{input_key}",
            )
    else:
        if st.button(load_label, key=f"load_example_{input_key}"):
            st.session_state[input_key] = example_text

    if help_text:
        st.caption(help_text)


def show_tip_box():
    st.sidebar.markdown("---")
    with st.sidebar.expander("☕ 打赏支持 / Buy me a coffee", expanded=False):
        st.markdown("如果这个工具帮到你，可以扫码支持一下开发维护。")
        if os.path.exists(TIP_CODE_PATH):
            st.image(TIP_CODE_PATH, caption="感谢支持 ❤️", use_container_width=True)
        else:
            st.warning("未找到打赏二维码图片，请检查路径：")
            st.code(TIP_CODE_PATH)


@st.cache_data(show_spinner=False)
def get_go_mapping_paths():
    return {
        "term2gene": os.path.join(GO_MAPPING_DIR, "TERM2GENE_protein_coding.tsv"),
        "term2name": os.path.join(GO_MAPPING_DIR, "TERM2NAME_protein_coding.tsv"),
        "metadata": os.path.join(GO_MAPPING_DIR, "wheat_go_metadata.tsv"),
        "background": os.path.join(GO_MAPPING_DIR, "wheat_protein_coding_genes.tsv"),
    }


@st.cache_data(show_spinner=False)
def get_kegg_mapping_paths():
    return {
        "gene2ko": os.path.join(KEGG_MAPPING_DIR, "gene2ko_clean.tsv"),
        "ko2pathway": os.path.join(KEGG_MAPPING_DIR, "kegg_ko2pathway.tsv"),
        "pathway2name": os.path.join(KEGG_MAPPING_DIR, "kegg_pathway2name.tsv"),
    }


@st.cache_data(show_spinner=False)
def run_cached_go_enrichment(
    gene_ids,
    term2gene_path,
    term2name_path,
    metadata_path,
    background_path,
    min_size,
    max_size,
    padj_cutoff,
):
    return run_go_enrichment(
        gene_list=list(gene_ids),
        term2gene_path=term2gene_path,
        term2name_path=term2name_path,
        metadata_path=metadata_path,
        background_path=background_path,
        min_size=min_size,
        max_size=max_size,
        padj_cutoff=padj_cutoff,
    )


@st.cache_data(show_spinner=False)
def run_cached_kegg_enrichment(
    gene_ids,
    gene2ko_path,
    ko2pathway_path,
    pathway2name_path,
    min_size,
    max_size,
    pvalue_cutoff,
):
    return run_kegg_enrichment(
        gene_list=list(gene_ids),
        gene2ko_path=gene2ko_path,
        ko2pathway_path=ko2pathway_path,
        pathway2name_path=pathway2name_path,
        min_size=min_size,
        max_size=max_size,
        pvalue_cutoff=pvalue_cutoff,
    )
