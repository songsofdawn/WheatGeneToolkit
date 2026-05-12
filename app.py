import streamlit as st

from app_shared import show_tip_box
from utils.db_query import DatabaseUnavailableError

TOOL_LABELS = [
    "ReadMe",
    "基因功能注释及三代基因号转换",
    "基因cDNA & CDS & protein sequences下载",
    "中国春同源基因检索（自身同源 + Fielder）",
    "Fielder 基因 → 启动子序列",
    "中国春启动子抓取",
    "JASPAR Plants 启动子 motif 分析",
    "火山图分析",
    "GO富集分析",
    "KEGG富集分析",
]

def render_selected_page(page_name: str):
    """按需导入页面模块，避免冷启动加载所有重型依赖。"""
    if page_name == "ReadMe":
        from sections import readme_page
        readme_page.render()
    elif page_name == "基因功能注释及三代基因号转换":
        from sections import gene_info_page
        gene_info_page.render()
    elif page_name == "基因cDNA & CDS & protein sequences下载":
        from sections import sequences_page
        sequences_page.render()
    elif page_name == "中国春同源基因检索（自身同源 + Fielder）":
        from sections import homolog_page
        homolog_page.render()
    elif page_name == "Fielder 基因 → 启动子序列":
        from sections import fielder_promoter_page
        fielder_promoter_page.render()
    elif page_name == "中国春启动子抓取":
        from sections import cs_promoter_page
        cs_promoter_page.render()
    elif page_name == "JASPAR Plants 启动子 motif 分析":
        from sections import motif_page
        motif_page.render()
    elif page_name == "火山图分析":
        from sections import volcano_page
        volcano_page.render()
    elif page_name == "GO富集分析":
        from sections import go_page
        go_page.render()
    elif page_name == "KEGG富集分析":
        from sections import kegg_page
        kegg_page.render()
    else:
        st.error(f"未知页面：{page_name}")


def main():
    st.set_page_config(page_title="WheatGeneToolkit 小麦基因批量处理工具", layout="wide")
    st.title("WheatGeneToolkit 小麦基因批量处理工具")

    tool = st.sidebar.radio("选择功能，初次使用请阅读 ReadMe", TOOL_LABELS)
    show_tip_box()
    try:
        render_selected_page(tool)
    except DatabaseUnavailableError as exc:
        st.error(str(exc))
        st.info("如果是在 Streamlit Cloud 部署，请确认仓库中包含 data/db/manifest.json 以及 data/db/ 下的 SQLite 分库文件。")


if __name__ == "__main__":
    main()
