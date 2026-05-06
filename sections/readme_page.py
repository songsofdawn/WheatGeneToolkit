import os

import streamlit as st

from app_shared import README_PATH


def render():
    st.header("说明文档 / ReadMe")

    if not os.path.exists(README_PATH):
        st.error("未找到 README.md")
        st.code(README_PATH)
        return

    with open(README_PATH, "r", encoding="utf-8") as f:
        st.markdown(f.read())

    st.info("请从左侧选择功能开始使用。")
