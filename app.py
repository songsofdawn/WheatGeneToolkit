import os
import pandas as pd
import streamlit as st

from utils.db_query import (
    get_primary_gene_id,
    get_gene_annotation,
    get_gene_core,
    get_sequences,
    get_promoter,
    get_fielder_best_hit,
    get_fielder_all_hits,
    get_fielder_promoter,
    get_cs_self_best_hit,
    get_cs_self_all_hits
)

from utils.go_enrichment import run_go_enrichment, create_go_barplot_bytes
from utils.kegg_enrichment import (
    run_kegg_enrichment,
    create_kegg_bubbleplot_bytes,
    create_kegg_barplot_bytes
)

# -------------------- 工具函数 --------------------
def read_gene_ids(uploaded_file, manual_input):
    if uploaded_file:
        return uploaded_file.read().decode("utf-8").splitlines()
    elif manual_input.strip():
        return manual_input.strip().splitlines()
    else:
        return []


# -------------------- 示例数据 --------------------
# 普通功能页的示例基因直接写在代码中，数量较少，主要用于演示输入格式。
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

# GO / KEGG 示例 DEG 从 data/ 目录中的外部 TXT 文件读取。
# 以后只需要修改这两个文件，不需要改 app.py：
#   data/example_go_deg_genes.txt
#   data/example_kegg_deg_genes.txt
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
GO_MAPPING_DIR = os.path.join(DATA_DIR, "go_mapping")
KEGG_MAPPING_DIR = os.path.join(DATA_DIR, "kegg_mapping")
TIP_CODE_PATH = os.path.join(DATA_DIR, "TipCode.jpg")
EXAMPLE_GO_DEG_PATH = os.path.join(DATA_DIR, "example_go_deg_genes.txt")
EXAMPLE_KEGG_DEG_PATH = os.path.join(DATA_DIR, "example_kegg_deg_genes.txt")

# 如果 data/ 目录中的示例文件不存在，使用下面的内置备用示例，避免网页报错。
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

def _as_text(items):
    return "\n".join(items)

def load_example_gene_text(file_path, fallback_items=None):
    """
    从外部 TXT 文件读取示例基因列表。

    文件格式：
    - 一行一个基因号；
    - 空行会自动忽略；
    - 以 # 开头的注释行会自动忽略。
    """
    if os.path.exists(file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            genes = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    continue
                genes.append(line)
        return _as_text(genes)

    if fallback_items is not None:
        return _as_text(fallback_items)

    return ""


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

def get_go_example_text():
    return load_example_gene_text(EXAMPLE_GO_DEG_PATH, DEFAULT_GO_DEG_GENES)

def get_kegg_example_text():
    return load_example_gene_text(EXAMPLE_KEGG_DEG_PATH, DEFAULT_KEGG_DEG_GENES)

def render_example_tools(
    input_key,
    example_text,
    load_label="🧪 加载示例数据",
    download_label=None,
    file_name=None,
    help_text=None,
):
    """
    在输入框上方显示示例数据按钮。

    - 点击“加载示例数据”后，会把示例内容写入对应 text_area 的 session_state；
    - 如果提供 download_label 和 file_name，则额外提供 TXT 示例文件下载。
    """
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
                key=f"download_example_{input_key}"
            )
    else:
        if st.button(load_label, key=f"load_example_{input_key}"):
            st.session_state[input_key] = example_text

    if help_text:
        st.caption(help_text)

# -------------------- 打赏提示 --------------------
def show_tip_box():
    """
    在侧边栏显示打赏二维码。
    图片路径：
    data/TipCode.jpg
    """
    st.sidebar.markdown("---")
    with st.sidebar.expander("☕ 打赏支持 / Buy me a coffee", expanded=False):
        st.markdown(
            """
            如果这个工具帮到了你，可以扫码支持一下开发维护。
            """
        )

        if os.path.exists(TIP_CODE_PATH):
            st.image(
                TIP_CODE_PATH,
                caption="感谢支持 ❤️",
                use_container_width=True
            )
        else:
            st.warning("未找到打赏二维码图片，请检查路径：")
            st.code(TIP_CODE_PATH)

# -------------------- 主程序 --------------------
def main():
    st.set_page_config(page_title="🌾 WheatGeneToolkit小麦基因批量处理工具", layout="wide")
    st.title("🌾 WheatGeneToolkit小麦基因批量处理工具")

    tool = st.sidebar.radio(
        "选择功能，初次使用请阅读ReadMe",
        [
            "ReadMe",
            "基因功能注释及三代基因号转换",
            "基因cDNA & CDS & protein sequences下载",
            "中国春同源基因检索（自身同源 + Fielder）",
            "Fielder 基因 → 启动子序列",
            "中国春启动子抓取",
            "GO富集分析",
            "KEGG富集分析"
        ]
    )
    show_tip_box()

    # -------------------- ReadMe --------------------
    if tool == "ReadMe":
        st.header("📖 WheatGeneToolkit / 小麦基因批量处理工具 - 使用说明")
        st.markdown("""
# WheatGeneToolkit / 小麦基因批量处理工具

🌾 **WheatGeneToolkit** is an integrated web-based toolkit for wheat gene annotation, gene ID conversion, sequence retrieval, homolog search, promoter extraction, and GO/KEGG enrichment analysis.  
🌾 **WheatGeneToolkit / 小麦基因批量处理工具** 是一个面向小麦功能基因组学研究的在线工具平台，集成了基因功能注释、基因号转换、序列下载、同源基因检索、启动子序列提取、GO 富集分析和 KEGG 富集分析等常用功能。

The current version is built with **Streamlit** and uses a **locally partitioned SQLite database** for stable and efficient online deployment.  
当前版本基于 **Streamlit** 构建，并采用 **本地分库 SQLite 数据库** 进行查询，避免直接加载大型数据库，更适合部署到 GitHub 和 Streamlit Cloud。

---

## Main Features / 主要功能

### 1. Gene annotation and gene ID conversion / 基因功能注释及基因号转换

Input Chinese Spring wheat gene IDs and retrieve basic annotation information.  
输入中国春小麦基因号后，可批量获取基础注释信息。

Supported outputs include / 支持输出内容包括：

- Input gene ID / 输入基因号
- Unified primary gene ID / 统一主键 `primary_gene_id`
- Third-generation gene ID / 三代基因号
- English functional description / 英文功能描述
- Chinese functional description / 中文功能描述

This module supports local alias-based gene ID conversion. If the input gene ID is an alias or older gene ID, the program will try to map it to the unified `primary_gene_id`.  
该模块支持基于本地 `gene_alias` 表的基因号转换。如果输入的是别名或旧版本基因号，程序会尝试自动转换为统一主键 `primary_gene_id`。

Example / 示例：

```text
TraesCS2D02G571200
TraesCS5B02G233300
TraesCS1A02G000100
```

---

### 2. cDNA / CDS / protein sequence retrieval / cDNA、CDS 和蛋白序列下载

Input Chinese Spring gene IDs and retrieve transcript-related sequences.  
输入中国春基因号后，可批量获取该基因对应转录本的序列信息。

Supported sequence types / 支持的序列类型：

- cDNA sequence / cDNA 序列
- CDS sequence / CDS 序列
- Protein sequence / 蛋白序列

The output is provided in standard FASTA format. One gene may correspond to multiple transcripts.  
输出结果为标准 FASTA 格式。同一个基因可能对应多个 transcript，因此可能返回多条序列。

---

### 3. Chinese Spring self-homolog search / 中国春自身同源基因检索

Input Chinese Spring gene IDs and retrieve homologous genes within the Chinese Spring genome.  
输入中国春基因号后，可查询中国春基因组内部的自身同源基因。

Supported outputs include / 支持输出内容包括：

- Best self-homolog hit / 最可信自身同源基因
- All candidate self-homolog genes / 全部候选自身同源基因
- Homolog type / 同源类型
- Confidence / 可信度
- Subgenome relationship / 亚基因组关系
- Chromosome relationship / 染色体关系
- Priority score / 优先级分数

This module is useful for analyzing gene duplication, homoeologous genes, and gene family relationships in hexaploid wheat.  
该功能适用于分析六倍体小麦中的基因复制、同源亚基因组基因和基因家族关系。

---

### 4. Chinese Spring to Fielder homolog search / 中国春到 Fielder 同源基因检索

Input Chinese Spring gene IDs and retrieve corresponding Fielder homologous genes.  
输入中国春基因号后，可查询对应的 Fielder 同源基因。

This function is especially useful for researchers working with Fielder as a transformation, genome editing, or functional validation background.  
该功能特别适合以 Fielder 为遗传转化、基因编辑或功能验证背景的研究者使用。

Supported outputs include / 支持输出内容包括：

- Best Fielder homolog / 最可信 Fielder 同源基因
- All candidate Fielder homologs / 全部候选 Fielder 同源基因
- Homolog type / 同源类型
- Confidence / 可信度
- Same subgenome information / 是否同亚基因组
- Same chromosome information / 是否同染色体
- Priority score / 优先级分数

---

### 5. Chinese Spring promoter extraction / 中国春启动子序列提取

Input Chinese Spring gene IDs and retrieve upstream promoter sequences.  
输入中国春基因号后，可批量获取其启动子序列。

Current promoter definition / 当前启动子定义：

```text
2000 bp upstream of ATG
ATG 上游 2000 bp
```

Features / 功能特点：

- Fixed promoter length of 2000 bp / 固定提取 ATG 上游 2000 bp
- Strand-aware extraction / 考虑正负链方向
- Reverse-complement handling for negative-strand genes / 对负链基因进行反向互补处理
- FASTA output / 支持 FASTA 格式下载
- Batch processing / 支持批量处理

For negative-strand genes, the output sequence is reverse-complemented and presented in the transcriptional 5' to 3' direction.  
对于负链基因，输出序列已经按照转录方向进行反向互补处理，最终结果统一为启动子 5' 到 3' 方向。也就是说，输出的启动子序列末端紧邻 ATG 所在位置。

---

### 6. Fielder promoter extraction / Fielder 启动子序列提取

Input Fielder gene IDs and retrieve Fielder promoter sequences.  
输入 Fielder 基因号后，可批量获取 Fielder 启动子序列。

Fielder gene IDs usually follow this format / Fielder 基因号通常为如下格式：

```text
TraesFLD5B01G105200
```

Current promoter definition / 当前启动子定义：

```text
2000 bp upstream of ATG
ATG 上游 2000 bp
```

Features / 功能特点：

- Supports `TraesFLD...` gene IDs / 支持 `TraesFLD...` 格式基因号
- Fixed promoter length of 2000 bp / 固定提取 ATG 上游 2000 bp
- Strand-aware extraction / 考虑正负链方向
- Reverse-complement handling for negative-strand genes / 对负链基因进行反向互补处理
- FASTA output / 支持 FASTA 格式下载

---

### 7. GO enrichment analysis / GO 富集分析

Input a DEG list with one gene ID per line and perform GO enrichment analysis.  
输入差异基因列表，一行一个基因号，即可进行 GO 富集分析。

The GO enrichment module outputs / GO 富集分析模块输出内容包括：

- Full GO enrichment result table / GO 富集结果总表
- Significant GO term table / 显著 GO 条目表
- GO enrichment bar plot / GO 富集条形图
- Analysis summary table / 分析摘要表

The enrichment analysis is based on local wheat GO annotation files and uses hypergeometric testing with multiple-testing correction.  
GO 富集分析基于本地小麦 GO 注释文件，使用超几何检验进行富集分析，并进行多重检验校正。

Local GO annotation files are stored in / 本地 GO 注释文件位于：

```text
data/go_mapping/
├── TERM2GENE_protein_coding.tsv
├── TERM2NAME_protein_coding.tsv
├── wheat_go_metadata.tsv
└── wheat_protein_coding_genes.tsv
```

GO enrichment workflow / GO 富集分析流程：

```text
Input DEG list
        ↓
Map genes to GO terms
        ↓
Filter GO terms by gene set size
        ↓
Hypergeometric enrichment test
        ↓
Multiple-testing correction
        ↓
Output result tables and bar plot
```

```text
输入差异基因列表
        ↓
将基因映射到 GO 条目
        ↓
按照基因集大小过滤 GO 条目
        ↓
进行超几何检验
        ↓
进行多重检验校正
        ↓
输出富集结果表和条形图
```

Notes / 注意事项：

- GO enrichment results depend on the completeness of the local GO annotation.  
  GO 富集结果依赖本地 GO 注释文件的完整性。
- Genes without GO annotation will not contribute to GO enrichment testing.  
  没有 GO 注释的基因不会进入 GO 富集检验。
- Enrichment results indicate over-representation of GO terms, not direct biological causality.  
  富集结果只能说明某些 GO 条目在输入基因集中显著偏多，不能直接证明因果关系。

---

### 8. KEGG enrichment analysis / KEGG 富集分析

Input a DEG list with one gene ID per line and perform KEGG pathway enrichment analysis.  
输入差异基因列表，一行一个基因号，即可进行 KEGG 通路富集分析。

The KEGG enrichment module outputs / KEGG 富集分析模块输出内容包括：

- Full KEGG enrichment result table / KEGG 富集结果总表
- Significant pathway table / 显著通路结果表
- KEGG bubble plot / KEGG 富集气泡图
- KEGG bar plot / KEGG 富集条形图
- Analysis summary table / KEGG 分析摘要表

This module is based on local gene-KO and KO-pathway mapping files. It does not depend on online KEGG queries during analysis.  
该模块基于本地 gene-KO 和 KO-pathway 映射文件进行分析，运行时不依赖在线 KEGG 查询。

Local KEGG files are stored in / 本地 KEGG 注释文件位于：

```text
data/kegg_mapping/
├── gene2ko_clean.tsv
├── kegg_ko2pathway.tsv
└── kegg_pathway2name.tsv
```

File descriptions / 文件说明：

**`gene2ko_clean.tsv`**  
Mapping between wheat genes and KEGG Orthology identifiers.  
小麦基因与 KEGG Orthology，也就是 KO 编号之间的映射表。

Example / 示例：

```text
TraesCS1A02G000100    K00001
TraesCS1A02G000200    K14488
```

**`kegg_ko2pathway.tsv`**  
Mapping between KO identifiers and KEGG pathway IDs.  
KO 编号与 KEGG pathway 编号之间的映射表。

Example / 示例：

```text
K00001    map00010
K14488    map04075
```

**`kegg_pathway2name.tsv`**  
Mapping between KEGG pathway IDs and pathway names.  
KEGG pathway 编号与通路名称之间的映射表。

Example / 示例：

```text
map00010    Glycolysis / Gluconeogenesis
map04075    Plant hormone signal transduction
```

KEGG enrichment workflow / KEGG 富集分析流程：

```text
Input DEG list
        ↓
Map wheat genes to KO identifiers
        ↓
Map KO identifiers to KEGG pathways
        ↓
Count pathway-level KO hits
        ↓
Perform hypergeometric enrichment test
        ↓
Output enrichment tables, bubble plot, and bar plot
```

```text
输入差异基因列表
        ↓
将小麦基因映射到 KO 编号
        ↓
将 KO 编号映射到 KEGG pathway
        ↓
统计每个 pathway 中被命中的 KO 数量
        ↓
进行超几何富集检验
        ↓
输出富集结果表、气泡图和条形图
```

Common output fields / 常见输出字段说明：

| Field | English description | 中文说明 |
|---|---|---|
| pathway_id | KEGG pathway ID | KEGG 通路编号 |
| pathway_name | KEGG pathway name | KEGG 通路名称 |
| k | Number of input KOs mapped to this pathway | 输入基因对应 KO 中命中该通路的 KO 数 |
| K | Number of background KOs in this pathway | 背景 KO 中属于该通路的 KO 数 |
| n | Number of input KOs used for enrichment | 输入基因中成功映射到 KO 的数量 |
| N | Number of background KOs used for enrichment | 背景中可用于分析的 KO 总数 |
| pvalue | P-value from hypergeometric test | 超几何检验得到的 P 值 |
| gene_ratio | Ratio of input KOs mapped to this pathway | 输入 KO 中命中该通路的比例 |
| background_ratio | Ratio of background KOs in this pathway | 背景 KO 中属于该通路的比例 |
| hit_ko | KO identifiers mapped to this pathway | 命中该通路的 KO 编号 |
| hit_genes | Input genes mapped to this pathway | 命中该通路的输入基因 |

Notes / 注意事项：

- KEGG enrichment is performed at the KO level rather than directly at the raw gene ID level.  
  KEGG 富集是在 KO 层面进行统计，而不是直接按原始基因号统计。
- Genes without KO annotation will be excluded from KEGG enrichment testing.  
  没有 KO 注释的基因不会进入 KEGG 富集检验。
- One gene may correspond to one or more KO identifiers.  
  一个基因可能对应一个或多个 KO 编号。
- One KO identifier may participate in multiple KEGG pathways.  
  一个 KO 编号也可能参与多个 KEGG pathway。
- KEGG enrichment indicates over-representation of pathways among input genes.  
  KEGG 富集结果说明某些通路在输入基因对应的 KO 中显著偏多。
- KEGG enrichment does not directly indicate whether a pathway is activated or repressed.  
  KEGG 富集结果不能直接说明通路被激活或被抑制。
- To infer pathway activation or repression, RNA-seq fold-change direction, expression pattern, and biological context should be considered together.  
  如果需要判断通路上调或下调，需要结合 RNA-seq 的 log2FoldChange、表达趋势和具体生物学背景进一步解释。

---

## Database Design / 数据库设计

The original large SQLite database has been partitioned into multiple smaller query-ready SQLite databases.  
原始大型 SQLite 数据库已经被拆分为多个可以直接查询的小型 SQLite 数据库。

This design avoids loading, merging, or decompressing a large database during web app startup.  
这种设计避免了网站启动时加载、合并或解压大型数据库，从而更适合 GitHub 和 Streamlit Cloud 部署。

Current database structure / 当前数据库结构：

```text
data/db/
├── manifest.json
├── core/
│   ├── fielder_gene_core.db
│   ├── gene_alias.db
│   ├── gene_annotation.db
│   ├── gene_core.db
│   └── transcript_core.db
│
├── gene_sequence_resource/
│   ├── gene_sequence_resource_1A.db
│   ├── gene_sequence_resource_1B.db
│   ├── gene_sequence_resource_1D.db
│   ├── ...
│   └── gene_sequence_resource_unknown.db
│
├── gene_promoter_sequence/
│   ├── gene_promoter_sequence_1A.db
│   ├── gene_promoter_sequence_1B.db
│   ├── gene_promoter_sequence_1D.db
│   ├── ...
│   └── gene_promoter_sequence_unknown.db
│
├── fielder_promoter_sequence/
│   ├── fielder_promoter_sequence_1A.db
│   ├── fielder_promoter_sequence_1B.db
│   ├── fielder_promoter_sequence_1D.db
│   ├── ...
│   └── fielder_promoter_sequence_unknown.db
│
├── gene_structure_feature/
│   ├── gene_structure_feature_1A.db
│   ├── gene_structure_feature_1B.db
│   ├── gene_structure_feature_1D.db
│   ├── ...
│   └── gene_structure_feature_unknown.db
│
└── homolog/
    ├── cs_self_homolog_map.db
    └── homolog_map.db
```

The `manifest.json` file records how each table is stored and how each SQLite shard should be located.  
`manifest.json` 文件记录了每张表的存储方式以及每个 SQLite 分库的位置。

During query execution, the program automatically identifies the chromosome from the gene ID.  
查询时，程序会自动从基因号中识别染色体。

Example / 示例：

```text
TraesCS5B02G233300 → 5B
TraesFLD5B01G105200 → 5B
```

Then the program opens only the corresponding small SQLite database.  
然后程序只打开对应染色体的小型 SQLite 数据库进行查询。

---

## Input Suggestions / 输入建议

- Use one gene ID per line.  
  建议每行输入一个基因号。
- Chinese Spring genes usually follow the `TraesCS...` format.  
  中国春基因通常使用 `TraesCS...` 格式。
- Fielder genes usually follow the `TraesFLD...` format.  
  Fielder 基因通常使用 `TraesFLD...` 格式。
- For batch analysis, uploading a TXT file is recommended.  
  批量分析时，推荐上传 TXT 文件。
- Transcript suffixes such as `.1` and `.2` are generally not required.  
  一般不需要添加 `.1`、`.2` 等 transcript 后缀。

---

## Notes / 注意事项

- This tool is designed for wheat gene list processing and functional genomics analysis.  
  本工具主要用于小麦基因列表处理和功能基因组学分析。
- The current promoter definition is fixed as 2000 bp upstream of ATG.  
  当前启动子定义固定为 ATG 上游 2000 bp。
- Homolog relationships are computationally inferred and should be interpreted with caution.  
  同源关系为计算推断结果，应谨慎解释。
- GO enrichment results depend on the quality and completeness of local GO annotation files.  
  GO 富集结果依赖本地 GO 注释文件的质量和完整性。
- KEGG enrichment results depend on the quality and completeness of local gene-KO and KO-pathway mapping files.  
  KEGG 富集结果依赖本地 gene-KO 和 KO-pathway 映射文件的质量和完整性。
- Enrichment results should be interpreted together with biological knowledge and experimental validation.  
  富集分析结果应结合生物学知识和实验验证进行综合判断。

---

## Developer / 开发人员

NWAFU 科研楼 2303 wyz  
GitHub: songsofdawn

---

## Disclaimer / 免责声明

WheatGeneToolkit is provided for research use only.  
WheatGeneToolkit 仅供科研使用。

The results generated by this platform should be interpreted together with original genome annotations, literature evidence, experimental data, and biological context.  
本平台生成的结果应结合原始基因组注释、文献证据、实验数据和具体生物学背景进行综合判断。
        """)
        st.info("📌 请从左侧选择功能开始使用 / Please select a function from the sidebar to get started.")
        return

    # -------------------- 功能 1 --------------------
    if tool == "基因功能注释及三代基因号转换":
        st.header("🔍 基因功能注释及三代基因号转换（本地数据库版）")

        uploaded_file = st.file_uploader(
            "上传 TXT 文件（基因号一行一个，不要加“.1”）",
            type=["txt"],
            key="file_gene_info"
        )
        render_example_tools(
            input_key="input_gene_info",
            example_text=_as_text(EXAMPLE_CS_GENES),
            load_label="🧪 加载示例中国春基因",
            help_text="示例会自动填入下方输入框，可直接点击开始查询。"
        )
        manual_input = st.text_area(
            "或者手动输入基因号（每行一个）",
            key="input_gene_info"
        )

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        st.info(f"待查询基因数：{len(gene_ids)}")

        if st.button("开始查询", key="btn_gene_info"):
            progress_bar = st.progress(0)
            status_text = st.empty()

            results = []

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在查询：{input_gene_id} ({idx}/{len(gene_ids)})")

                primary_gene_id = get_primary_gene_id(input_gene_id)
                anno_df = get_gene_annotation(input_gene_id)
                core_df = get_gene_core(input_gene_id)

                third_id = None
                description_en = None
                description_zh = None

                if not anno_df.empty:
                    description_en = anno_df.iloc[0].get("description_en", None)
                    description_zh = anno_df.iloc[0].get("description_zh", None)

                if not core_df.empty:
                    third_id = core_df.iloc[0].get("gene_id_v3", None)

                results.append([
                    input_gene_id,
                    primary_gene_id if primary_gene_id is not None else "未找到",
                    third_id if pd.notna(third_id) else "",
                    description_en if pd.notna(description_en) else "",
                    description_zh if pd.notna(description_zh) else ""
                ])

                progress_bar.progress(idx / len(gene_ids))

            df = pd.DataFrame(
                results,
                columns=["输入基因号", "统一主键", "三代基因号", "功能描述（英文）", "功能描述（中文）"]
            )

            st.success("✅ 查询完成！")
            st.dataframe(df, use_container_width=True)

            st.download_button(
                "📥 下载结果 CSV",
                df.to_csv(index=False).encode("utf-8-sig"),
                "gene_info_from_sqlite.csv",
                "text/csv"
            )

            status_text.empty()
            progress_bar.empty()

    # -------------------- 功能 2 --------------------
    elif tool == "基因cDNA & CDS & protein sequences下载":
        st.header("📍 cDNA / CDS / Protein 下载（本地数据库版）")

        uploaded_file = st.file_uploader("上传 TXT 文件（基因号）", type=["txt"], key="file_sequences")
        render_example_tools(
            input_key="input_sequences",
            example_text=_as_text(EXAMPLE_CS_GENES),
            load_label="🧪 加载示例中国春基因",
            help_text="示例会自动填入下方输入框，可直接获取 cDNA / CDS / protein。"
        )
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_sequences")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

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
                    len(protein_df)
                ])

                for _, row in cdna_df.iterrows():
                    tx_id = row.get("transcript_id", "")
                    seq = row.get("sequence", "")
                    header = f">{input_gene_id}|{primary_gene_id}|{tx_id}|cdna"
                    cdna_records.append(f"{header}\n{seq}\n")

                for _, row in cds_df.iterrows():
                    tx_id = row.get("transcript_id", "")
                    seq = row.get("sequence", "")
                    header = f">{input_gene_id}|{primary_gene_id}|{tx_id}|cds"
                    cds_records.append(f"{header}\n{seq}\n")

                for _, row in protein_df.iterrows():
                    tx_id = row.get("transcript_id", "")
                    seq = row.get("sequence", "")
                    header = f">{input_gene_id}|{primary_gene_id}|{tx_id}|protein"
                    protein_records.append(f"{header}\n{seq}\n")

                if len(cdna_df) == 0 and len(cds_df) == 0 and len(protein_df) == 0:
                    failed_genes.append(input_gene_id)

                progress.progress(idx / len(gene_ids))

            summary_df = pd.DataFrame(
                summary_rows,
                columns=["输入基因号", "统一主键", "cDNA条数", "CDS条数", "Protein条数"]
            )

            st.success("✅ 序列查询完成！")
            st.dataframe(summary_df, use_container_width=True)

            st.download_button("📥 下载 cDNA FASTA", "".join(cdna_records), "cdna_sequences.fasta", "text/plain")
            st.download_button("📥 下载 CDS FASTA", "".join(cds_records), "cds_sequences.fasta", "text/plain")
            st.download_button("📥 下载 Protein FASTA ", "".join(protein_records), "protein_sequences.fasta", "text/plain")
            st.download_button(
                "📥 下载序列统计 CSV",
                summary_df.to_csv(index=False).encode("utf-8-sig"),
                "sequence_summary.csv",
                "text/csv"
            )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到序列: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 3 --------------------
    elif tool == "中国春同源基因检索（自身同源 + Fielder）":
        st.header("🧬 中国春同源基因检索（自身同源 + Fielder）")

        uploaded_file = st.file_uploader("上传 TXT 文件（基因号一行一个）", type=["txt"], key="file_blast")
        render_example_tools(
            input_key="input_blast",
            example_text=_as_text(EXAMPLE_CS_GENES),
            load_label="🧪 加载示例中国春基因",
            help_text="示例会自动填入下方输入框，可直接测试自身同源和 Fielder 同源查询。"
        )
        manual_input = st.text_area("或者手动输入中国春基因号（每行一个）", height=200, key="input_blast")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("💡 请上传文件或输入中国春基因号")
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
                        r.get("priority_score", "")
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
                        r.get("priority_score", "")
                    ])
                    tmp = fld_all_df_one.copy()
                    tmp.insert(0, "input_gene_id", input_gene_id)
                    fld_all_rows.append(tmp)

                progress.progress(idx / len(gene_ids))

            self_best_result_df = pd.DataFrame(
                self_best_rows,
                columns=["输入基因号", "中国春主键", "默认中国春自身同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"]
            )

            fld_best_result_df = pd.DataFrame(
                fld_best_rows,
                columns=["输入基因号", "中国春主键", "默认Fielder同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"]
            )

            self_all_df = pd.concat(self_all_rows, ignore_index=True) if self_all_rows else pd.DataFrame()
            fld_all_df = pd.concat(fld_all_rows, ignore_index=True) if fld_all_rows else pd.DataFrame()

            st.success("✅ 查询完成！")

            st.subheader("中国春自身同源：默认最可信结果")
            st.dataframe(self_best_result_df, use_container_width=True)

            st.subheader("中国春自身同源：全部候选")
            if not self_all_df.empty:
                st.dataframe(self_all_df, use_container_width=True)
                st.download_button(
                    "📥 下载中国春自身同源全部候选 CSV",
                    self_all_df.to_csv(index=False).encode("utf-8-sig"),
                    "cs_self_homolog_all_hits.csv",
                    "text/csv"
                )

            st.subheader("Fielder 同源：默认最可信结果")
            st.dataframe(fld_best_result_df, use_container_width=True)

            st.subheader("Fielder 同源：全部候选")
            if not fld_all_df.empty:
                st.dataframe(fld_all_df, use_container_width=True)
                st.download_button(
                    "📥 下载 Fielder 同源全部候选 CSV",
                    fld_all_df.to_csv(index=False).encode("utf-8-sig"),
                    "cs_to_fielder_all_hits.csv",
                    "text/csv"
                )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未找到同源结果: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 4 --------------------
    elif tool == "Fielder 基因 → 启动子序列":
        st.header("🌱 Fielder 基因 → 启动子序列（本地数据库版）")

        uploaded_file = st.file_uploader("上传 Fielder 基因列表 TXT", type=["txt"], key="file_fielder_promoter")
        render_example_tools(
            input_key="input_fielder_promoter",
            example_text=_as_text(EXAMPLE_FIELDER_GENES),
            load_label="🧪 加载示例 Fielder 基因",
            help_text="示例会自动填入下方输入框，用于测试 Fielder 启动子提取。"
        )
        manual_input = st.text_area("或者手动输入 Fielder 基因号（每行一个）", key="input_fielder_promoter")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入 Fielder 基因号")
            st.stop()

        st.info("当前数据库内置的是 promoter_2000，因此这里固定返回 ATG 上游 2000 bp 启动子序列。")
        st.info(f"待处理基因数：{len(gene_ids)}")

        if st.button("开始抓取 Fielder 启动子", key="btn_fielder_promoter"):
            progress = st.progress(0)
            status_text = st.empty()

            fasta_records = []
            summary_rows = []
            failed_genes = []

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在查询: {input_gene_id} ({idx}/{len(gene_ids)})")

                promoter_df = get_fielder_promoter(input_gene_id)

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

                header = f">{input_gene_id}|{primary_id}|promoter_2000|{chrom}:{promoter_start}-{promoter_end}|{strand}"
                fasta_records.append(f"{header}\n{promoter_seq}\n")

                summary_rows.append([
                    input_gene_id,
                    primary_id if pd.notna(primary_id) else "",
                    chrom if pd.notna(chrom) else "",
                    promoter_start if pd.notna(promoter_start) else "",
                    promoter_end if pd.notna(promoter_end) else "",
                    strand if pd.notna(strand) else "",
                    promoter_length if pd.notna(promoter_length) else 0
                ])

                progress.progress(idx / len(gene_ids))

            summary_df = pd.DataFrame(
                summary_rows,
                columns=["输入基因号", "Fielder主键", "染色体", "启动子起点", "启动子终点", "链方向", "启动子长度"]
            )

            st.success("🎉 Fielder 启动子抓取完成！")
            st.dataframe(summary_df, use_container_width=True)

            st.download_button(
                "📥 下载 Fielder 启动子 FASTA",
                "".join(fasta_records),
                "fielder_promoter_sequences.fasta",
                "text/plain"
            )

            st.download_button(
                "📥 下载 Fielder 启动子统计 CSV",
                summary_df.to_csv(index=False).encode("utf-8-sig"),
                "fielder_promoter_summary.csv",
                "text/csv"
            )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到启动子序列: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 5 --------------------
    elif tool == "中国春启动子抓取":
        st.header("🌱 中国春基因启动子抓取（本地数据库版）")

        uploaded_file = st.file_uploader("上传基因列表 TXT", type=["txt"], key="file_cs_promoter")
        render_example_tools(
            input_key="input_cs_promoter",
            example_text=_as_text(EXAMPLE_CS_GENES),
            load_label="🧪 加载示例中国春基因",
            help_text="示例会自动填入下方输入框，用于测试中国春启动子提取。"
        )
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_cs_promoter")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        st.info("当前数据库内置的是 promoter_2000，因此这里固定返回 ATG 上游 2000 bp 启动子序列。")
        st.info(f"待处理基因数: {len(gene_ids)}")

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

                header = f">{input_gene_id}|{primary_id}|promoter_2000|{chrom}:{promoter_start}-{promoter_end}|{strand}"
                fasta_records.append(f"{header}\n{promoter_seq}\n")

                summary_rows.append([
                    input_gene_id,
                    primary_id if pd.notna(primary_id) else "",
                    chrom if pd.notna(chrom) else "",
                    promoter_start if pd.notna(promoter_start) else "",
                    promoter_end if pd.notna(promoter_end) else "",
                    strand if pd.notna(strand) else "",
                    promoter_length if pd.notna(promoter_length) else 0
                ])

                progress.progress(idx / len(gene_ids))

            summary_df = pd.DataFrame(
                summary_rows,
                columns=["输入基因号", "统一主键", "染色体", "启动子起点", "启动子终点", "链方向", "启动子长度"]
            )

            st.success("🎉 启动子抓取完成！")
            st.dataframe(summary_df, use_container_width=True)

            st.download_button(
                "📥 下载启动子 FASTA",
                "".join(fasta_records),
                "cs_promoter_sequences.fasta",
                "text/plain"
            )

            st.download_button(
                "📥 下载启动子统计 CSV",
                summary_df.to_csv(index=False).encode("utf-8-sig"),
                "cs_promoter_summary.csv",
                "text/csv"
            )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到启动子序列: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 6 --------------------
    elif tool == "GO富集分析":
        st.header("📊 GO 富集分析")
        st.caption("输入 DEG 列表（一行一个基因号），输出 GO 富集结果表和条形图。")

        uploaded_file = st.file_uploader("上传 DEG TXT 文件（一行一个基因号）", type=["txt"], key="file_go_deg")
        render_example_tools(
            input_key="input_go_deg",
            example_text=get_go_example_text(),
            load_label="🧪 加载 GO 示例 DEG",
            download_label="📄 下载 GO 示例 DEG 文件",
            file_name="example_go_deg_genes.txt",
            help_text="示例 DEG 从 data/example_go_deg_genes.txt 读取；你可以直接修改该文件来更新 GO 示例。"
        )
        manual_input = st.text_area("或者手动输入 DEG（每行一个）", height=200, key="input_go_deg")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传 DEG 文件或手动输入基因号")
            st.stop()

        col1, col2, col3 = st.columns(3)
        with col1:
            top_n = st.number_input("每个大类展示前 N 个 term", min_value=3, max_value=30, value=10, step=1)
        with col2:
            padj_cutoff = st.number_input("FDR 阈值", min_value=0.0001, max_value=1.0, value=0.05, step=0.01, format="%.4f")
        with col3:
            plot_metric = st.selectbox(
                "作图横轴",
                options=["fdr", "pvalue"],
                format_func=lambda x: "-log10(FDR)" if x == "fdr" else "-log10(P-value)"
            )

        min_size = st.number_input("最小 GO 基因集大小", min_value=1, max_value=50, value=3, step=1)
        max_size = st.number_input("最大 GO 基因集大小", min_value=10, max_value=10000, value=2000, step=10)

        go_paths = get_go_mapping_paths()
        term2gene_path = go_paths["term2gene"]
        term2name_path = go_paths["term2name"]
        metadata_path = go_paths["metadata"]
        background_path = go_paths["background"]

        required_files = [term2gene_path, term2name_path, metadata_path, background_path]
        missing_files = [fp for fp in required_files if not os.path.exists(fp)]

        if missing_files:
            st.error("以下 GO 注释文件不存在，请检查 data/go_mapping/ 目录：")
            for fp in missing_files:
                st.code(fp)
            st.stop()

        st.info(f"待分析基因数: {len(gene_ids)}")

        if st.button("开始 GO 富集分析", key="btn_go_enrichment"):
            with st.spinner("正在进行 GO 富集分析，请稍候..."):
                try:
                    results_df, sig_df, summary_df = run_cached_go_enrichment(
                        gene_ids=tuple(gene_ids),
                        term2gene_path=term2gene_path,
                        term2name_path=term2name_path,
                        metadata_path=metadata_path,
                        background_path=background_path,
                        min_size=min_size,
                        max_size=max_size,
                        padj_cutoff=padj_cutoff
                    )

                    if results_df.empty:
                        st.warning("没有得到任何 GO 富集结果，请检查输入基因 ID 是否与背景一致。")
                        st.stop()

                    st.success("✅ GO 富集分析完成！")

                    st.subheader("分析摘要")
                    st.dataframe(summary_df, use_container_width=True)

                    st.subheader("显著富集结果")
                    if sig_df.empty:
                        st.warning("当前 FDR 阈值下没有显著富集条目，下面展示全部结果。")
                        st.dataframe(results_df, use_container_width=True)
                    else:
                        st.dataframe(sig_df, use_container_width=True)

                    plot_df = sig_df if not sig_df.empty else results_df
                    plot_bytes = create_go_barplot_bytes(
                        df=plot_df,
                        top_n=top_n,
                        plot_metric=plot_metric
                    )

                    if plot_bytes is not None:
                        st.subheader("GO 富集条形图")
                        st.image(plot_bytes, caption="GO enrichment barplot", use_container_width=True)
                        st.download_button(
                            "📥 下载富集图 PNG",
                            data=plot_bytes,
                            file_name="GO_enrichment_barplot.png",
                            mime="image/png"
                        )

                    st.download_button(
                        "📥 下载全部结果 TSV",
                        data=results_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="GO_enrichment_results.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载显著结果 TSV",
                        data=sig_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="GO_enrichment_results_sig.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载分析摘要 TSV",
                        data=summary_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="GO_enrichment_summary.tsv",
                        mime="text/tab-separated-values"
                    )

                except Exception as e:
                    st.error(f"GO 富集分析失败：{e}")

    # -------------------- 功能 7 --------------------
    elif tool == "KEGG富集分析":
        st.header("📊 KEGG 富集分析")
        st.caption("输入 DEG 列表（一行一个基因号），基于本地 gene-KO 和 KEGG KO-pathway 映射进行 KEGG 富集分析。")

        uploaded_file = st.file_uploader(
            "上传 DEG TXT 文件（一行一个基因号）",
            type=["txt"],
            key="file_kegg_deg"
        )
        render_example_tools(
            input_key="input_kegg_deg",
            example_text=get_kegg_example_text(),
            load_label="🧪 加载 KEGG 示例 DEG",
            download_label="📄 下载 KEGG 示例 DEG 文件",
            file_name="example_kegg_deg_genes.txt",
            help_text="示例 DEG 从 data/example_kegg_deg_genes.txt 读取；你可以直接修改该文件来更新 KEGG 示例。"
        )

        manual_input = st.text_area(
            "或者手动输入 DEG（每行一个）",
            height=200,
            key="input_kegg_deg"
        )

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传 DEG 文件或手动输入基因号")
            st.stop()

        st.info(f"待分析基因数: {len(gene_ids)}")

        col1, col2, col3 = st.columns(3)

        with col1:
            top_n = st.number_input(
                "图中展示前 N 个通路",
                min_value=3,
                max_value=50,
                value=20,
                step=1
            )

        with col2:
            pvalue_cutoff = st.number_input(
                "P-value 阈值",
                min_value=0.0001,
                max_value=1.0,
                value=0.05,
                step=0.01,
                format="%.4f"
            )

        with col3:
            use_sig_only = st.checkbox(
                "优先展示显著通路",
                value=True,
                help="如果存在 pvalue 小于阈值的通路，则优先使用显著通路绘图；否则展示全部结果中 pvalue 最小的通路。"
            )

        col4, col5 = st.columns(2)

        with col4:
            min_size = st.number_input(
                "最小 pathway KO 数",
                min_value=1,
                max_value=100,
                value=3,
                step=1
            )

        with col5:
            max_size = st.number_input(
                "最大 pathway KO 数",
                min_value=10,
                max_value=5000,
                value=500,
                step=10
            )

        with st.expander("高级绘图参数"):
            col6, col7, col8 = st.columns(3)

            with col6:
                clip_minus_log10_p = st.checkbox(
                    "截断极端 -log10(P-value)",
                    value=True,
                    help="建议开启，避免极端小 pvalue 把横轴拉得过长。"
                )

            with col7:
                clip_quantile = st.number_input(
                    "截断分位数",
                    min_value=0.50,
                    max_value=1.00,
                    value=0.95,
                    step=0.01,
                    format="%.2f"
                )

            with col8:
                label_wrap_width = st.number_input(
                    "通路名称换行宽度",
                    min_value=20,
                    max_value=80,
                    value=42,
                    step=2
                )

            col9, col10 = st.columns(2)

            with col9:
                bubble_min_size = st.number_input(
                    "最小气泡大小",
                    min_value=5,
                    max_value=200,
                    value=25,
                    step=5
                )

            with col10:
                bubble_max_size = st.number_input(
                    "最大气泡大小",
                    min_value=50,
                    max_value=800,
                    value=220,
                    step=10
                )

        kegg_paths = get_kegg_mapping_paths()
        gene2ko_path = kegg_paths["gene2ko"]
        ko2pathway_path = kegg_paths["ko2pathway"]
        pathway2name_path = kegg_paths["pathway2name"]

        required_files = [
            gene2ko_path,
            ko2pathway_path,
            pathway2name_path
        ]

        missing_files = [fp for fp in required_files if not os.path.exists(fp)]

        if missing_files:
            st.error("以下 KEGG 注释文件不存在，请检查 data/kegg_mapping/ 目录：")
            for fp in missing_files:
                st.code(fp)
            st.stop()

        with st.expander("当前使用的 KEGG 本地注释文件"):
            st.code(gene2ko_path)
            st.code(ko2pathway_path)
            st.code(pathway2name_path)

        if st.button("开始 KEGG 富集分析", key="btn_kegg_enrichment"):
            with st.spinner("正在进行 KEGG 富集分析，请稍候..."):
                try:
                    results_df, sig_df, summary_df = run_cached_kegg_enrichment(
                        gene_ids=tuple(gene_ids),
                        gene2ko_path=gene2ko_path,
                        ko2pathway_path=ko2pathway_path,
                        pathway2name_path=pathway2name_path,
                        min_size=min_size,
                        max_size=max_size,
                        pvalue_cutoff=pvalue_cutoff
                    )

                    if results_df.empty:
                        st.warning("没有得到任何 KEGG 富集结果，请检查输入基因 ID 是否与 gene2ko_clean.tsv 中的基因 ID 一致。")
                        st.subheader("分析摘要")
                        st.dataframe(summary_df, use_container_width=True)
                        st.stop()

                    st.success("✅ KEGG 富集分析完成！")

                    st.subheader("分析摘要")
                    st.dataframe(summary_df, use_container_width=True)

                    st.subheader("显著富集结果")
                    if sig_df.empty:
                        st.warning("当前 P-value 阈值下没有显著富集通路，下面展示全部结果。")
                        st.dataframe(results_df, use_container_width=True)
                    else:
                        st.dataframe(sig_df, use_container_width=True)

                    if use_sig_only and not sig_df.empty:
                        plot_df = sig_df
                    else:
                        plot_df = results_df

                    bubble_bytes = create_kegg_bubbleplot_bytes(
                        df=plot_df,
                        top_n=top_n,
                        bubble_min_size=bubble_min_size,
                        bubble_max_size=bubble_max_size,
                        clip_minus_log10_p=clip_minus_log10_p,
                        clip_mode="quantile",
                        clip_quantile=clip_quantile,
                        clip_fixed_value=30,
                        min_x_cap=5,
                        label_wrap_width=label_wrap_width
                    )

                    barplot_bytes = create_kegg_barplot_bytes(
                        df=plot_df,
                        top_n=top_n,
                        clip_minus_log10_p=clip_minus_log10_p,
                        clip_mode="quantile",
                        clip_quantile=clip_quantile,
                        clip_fixed_value=30,
                        min_x_cap=5,
                        label_wrap_width=label_wrap_width
                    )

                    if bubble_bytes is not None:
                        st.subheader("KEGG 富集气泡图")
                        st.image(
                            bubble_bytes,
                            caption="KEGG enrichment bubble plot",
                            use_container_width=True
                        )

                        st.download_button(
                            "📥 下载 KEGG 气泡图 PNG",
                            data=bubble_bytes,
                            file_name="KEGG_enrichment_bubbleplot.png",
                            mime="image/png"
                        )

                    if barplot_bytes is not None:
                        st.subheader("KEGG 富集条形图")
                        st.image(
                            barplot_bytes,
                            caption="KEGG enrichment barplot",
                            use_container_width=True
                        )

                        st.download_button(
                            "📥 下载 KEGG 条形图 PNG",
                            data=barplot_bytes,
                            file_name="KEGG_enrichment_barplot.png",
                            mime="image/png"
                        )

                    st.download_button(
                        "📥 下载全部 KEGG 富集结果 TSV",
                        data=results_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="KEGG_enrichment_results_all.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载显著 KEGG 富集结果 TSV",
                        data=sig_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="KEGG_enrichment_results_sig.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载 KEGG 分析摘要 TSV",
                        data=summary_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="KEGG_enrichment_summary.tsv",
                        mime="text/tab-separated-values"
                    )

                except Exception as e:
                    st.error(f"KEGG 富集分析失败：{e}")

if __name__ == "__main__":
    main()
