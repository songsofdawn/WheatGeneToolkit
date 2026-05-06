# App Structure Guide

这个项目现在采用了“单入口 + 分页面模块”的结构：

```text
app.py
app_shared.py
sections/
  readme_page.py
  gene_info_page.py
  sequences_page.py
  homolog_page.py
  fielder_promoter_page.py
  cs_promoter_page.py
  go_page.py
  kegg_page.py
utils/
data/
```

## 各文件职责

- `app.py`
  只负责：
  1. 设置 Streamlit 页面参数
  2. 渲染侧边栏功能菜单
  3. 根据菜单分发到对应页面模块

- `app_shared.py`
  放公共内容：
  1. 路径常量
  2. 示例数据
  3. 通用 helper
  4. Streamlit 缓存包装
  5. 侧边栏打赏框

- `sections/*.py`
  每个文件只负责一个功能页的 `render()`。

- `utils/*.py`
  放底层数据查询和分析逻辑，不直接负责页面渲染。

## 以后新增功能的推荐方式

假设你要新增一个页面，比如“基因结构查询”。

### 1. 新建页面模块

在 `sections/` 下新增：

```text
sections/gene_structure_page.py
```

模块格式建议保持一致：

```python
import streamlit as st

def render():
    st.header("基因结构查询")
    # 这里写页面输入、按钮、结果展示
```

### 2. 如果有公共常量或示例数据，放到 `app_shared.py`

例如：

- 新示例基因列表
- 新数据文件路径
- 新缓存函数

都优先放到 `app_shared.py`，不要散落在多个页面文件里。

### 3. 如果有底层查询逻辑，放到 `utils/`

原则是：

- 页面交互写在 `sections/`
- 数据处理和查询写在 `utils/`

这样后面排查问题时会更清楚。

### 4. 在 `app.py` 注册新页面

你只需要改两处：

1. `TOOL_LABELS` 增加一个菜单名
2. `PAGE_RENDERERS` 增加一个映射

示例：

```python
from sections import gene_structure_page
```

然后：

```python
TOOL_LABELS = [
    ...,
    "基因结构查询",
]
```

```python
PAGE_RENDERERS = {
    ...,
    "基因结构查询": gene_structure_page.render,
}
```

## 当前这次拆分带来的直接收益

- `app.py` 不再承载所有页面逻辑，主入口更清晰
- 每个功能页都可以单独维护
- 新功能只需要加新页面文件，不需要继续把 `app.py` 堆大
- 更容易定位 bug，修改范围更小
- 后续如果要继续拆分页面内部逻辑，也更容易做

## 后续建议

如果你继续扩展这个项目，建议按下面的方向继续保持：

1. 页面文件只保留 UI 和交互流程
2. 分析逻辑尽量下沉到 `utils/`
3. 复用的路径、示例、缓存统一放 `app_shared.py`
4. 大段文档内容尽量放到独立文件中，而不是直接塞进页面代码

这样项目会比较稳，也更适合后续部署和持续维护。
