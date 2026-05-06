import streamlit as st

from app_shared import show_tip_box
from sections import (
    cs_promoter_page,
    fielder_promoter_page,
    gene_info_page,
    go_page,
    homolog_page,
    kegg_page,
    readme_page,
    sequences_page,
)

TOOL_LABELS = [
    "ReadMe",
    "基因功能注释及三代基因号转换",
    "基因cDNA & CDS & protein sequences下载",
    "中国春同源基因检索（自身同源 + Fielder）",
    "Fielder 基因 → 启动子序列",
    "中国春启动子抓取",
    "GO富集分析",
    "KEGG富集分析",
]

PAGE_RENDERERS = {
    "ReadMe": readme_page.render,
    "基因功能注释及三代基因号转换": gene_info_page.render,
    "基因cDNA & CDS & protein sequences下载": sequences_page.render,
    "中国春同源基因检索（自身同源 + Fielder）": homolog_page.render,
    "Fielder 基因 → 启动子序列": fielder_promoter_page.render,
    "中国春启动子抓取": cs_promoter_page.render,
    "GO富集分析": go_page.render,
    "KEGG富集分析": kegg_page.render,
}


def main():
    st.set_page_config(page_title="WheatGeneToolkit 小麦基因批量处理工具", layout="wide")
    st.title("WheatGeneToolkit 小麦基因批量处理工具")

    tool = st.sidebar.radio("选择功能，初次使用请阅读 ReadMe", TOOL_LABELS)
    show_tip_box()
    PAGE_RENDERERS[tool]()


if __name__ == "__main__":
    main()
