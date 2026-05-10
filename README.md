# WheatGeneToolkit / 小麦基因批量处理工具

🌾 **WheatGeneToolkit** is an integrated web-based toolkit for wheat gene annotation, gene ID conversion, sequence retrieval, homolog search, promoter extraction, JASPAR Plants promoter motif scanning, and GO/KEGG enrichment analysis.

🌾 **WheatGeneToolkit / 小麦基因批量处理工具** 是一个面向小麦功能基因组学研究的在线工具平台，集成了基因功能注释、基因号转换、序列下载、同源基因检索、启动子序列提取、JASPAR Plants 启动子 motif 扫描、GO 富集分析和 KEGG 富集分析等常用功能。

This platform is designed for wheat researchers working on RNA-seq, GWAS, QTL mapping, gene family analysis, promoter analysis, and molecular biology experiments.

本工具适用于 RNA-seq、GWAS、QTL 定位、基因家族分析、启动子分析和分子生物学实验等小麦研究场景，尤其适合对大批量基因列表进行快速处理。

The current version is built with **Streamlit** and uses a **locally partitioned SQLite database** for stable and efficient online deployment.

当前版本基于 **Streamlit** 构建，并采用 **本地分库 SQLite 数据库** 进行查询，避免直接加载大型数据库，更适合部署到 GitHub 和 Streamlit Cloud。

Pre-launch notes:

- The production query backend is `data/db/manifest.json` plus the partitioned SQLite files under `data/db/`; the old `wheat_toolkit.db` is only a source/legacy construction file.
- If the database bundle is missing, the app will show a friendly message asking you to check `data/db/manifest.json` and the SQLite shard files.
- Batch gene input is supported for gene info, sequence, homolog, promoter, GO, and KEGG modules. Very large inputs may take longer; test with a small list before formal analysis.
- Example Chinese Spring IDs: `TraesCS2D02G571200`, `TraesCS5B02G233300`, `TraesCS6A02G189300`. Example Fielder ID: `TraesFLD5B01G105200`.
- Empty sequence/promoter results do not create FASTA downloads; check whether the gene ID exists, whether the species/version matches, and whether the database bundle is complete.

上线前提示：

- 当前线上查询后端是 `data/db/manifest.json` 和 `data/db/` 下的 SQLite 分库；旧的 `wheat_toolkit.db` 仅作为历史单库或构建来源。
- 如果数据库文件缺失，网站会提示检查 `data/db/manifest.json` 和 SQLite 分库文件。
- gene info、序列、同源、启动子、GO、KEGG 均支持批量输入；输入很多基因时可能需要等待，正式分析前建议先用少量基因测试。
- 示例中国春基因：`TraesCS2D02G571200`、`TraesCS5B02G233300`、`TraesCS6A02G189300`。示例 Fielder 基因：`TraesFLD5B01G105200`。
- 序列或启动子查询没有成功结果时不会生成 FASTA 下载文件，请检查基因 ID、物种版本和数据库文件是否完整。

---

## Online App / 在线访问

The app can be deployed on Streamlit Cloud.

本项目可部署到 Streamlit Cloud 进行在线访问。

```text
Main file path: app.py
Python version: 3.11
```

---

## Main Features / 主要功能

### 1. Gene annotation and gene ID conversion / 基因功能注释及基因号转换

Input Chinese Spring wheat gene IDs and retrieve basic annotation information.

输入中国春小麦基因号后，可批量获取基础注释信息。

Supported outputs include:

支持输出内容包括：

- Input gene ID / 输入基因号
- Unified primary gene ID / 统一主键 `primary_gene_id`
- Third-generation gene ID / 三代基因号
- English functional description / 英文功能描述
- Chinese functional description / 中文功能描述

This module supports local alias-based gene ID conversion. If the input gene ID is an alias or older gene ID, the program will try to map it to the unified `primary_gene_id`.

该模块支持基于本地 `gene_alias` 表的基因号转换。如果输入的是别名或旧版本基因号，程序会尝试自动转换为统一主键 `primary_gene_id`。

Example input:

示例输入：

```text
TraesCS2D02G571200
TraesCS5B02G233300
TraesCS1A02G000100
```

---

### 2. cDNA / CDS / protein sequence retrieval / cDNA、CDS 和蛋白序列下载

Input Chinese Spring gene IDs and retrieve transcript-related sequences.

输入中国春基因号后，可批量获取该基因对应转录本的序列信息。

Supported sequence types:

支持的序列类型包括：

- cDNA sequence / cDNA 序列
- CDS sequence / CDS 序列
- Protein sequence / 蛋白序列

The output is provided in standard FASTA format. One gene may correspond to multiple transcripts.

输出结果为标准 FASTA 格式。同一个基因可能对应多个 transcript，因此可能返回多条序列。

---

### 3. Chinese Spring self-homolog search / 中国春自身同源基因检索

Input Chinese Spring gene IDs and retrieve homologous genes within the Chinese Spring genome.

输入中国春基因号后，可查询中国春基因组内部的自身同源基因。

Supported outputs include:

支持输出内容包括：

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

Supported outputs include:

支持输出内容包括：

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

Current promoter definition:

当前启动子定义：

```text
2000 bp upstream of ATG
```

Features:

功能特点：

- Fixed promoter length of 2000 bp / 固定提取 ATG 上游 2000 bp
- Strand-aware extraction / 考虑正负链方向
- Reverse-complement handling for negative-strand genes / 对负链基因进行反向互补处理
- FASTA output / 支持 FASTA 格式下载
- Batch processing / 支持批量处理

For negative-strand genes, the output sequence is reverse-complemented and presented in the transcriptional 5' to 3' direction.

对于负链基因，输出序列已经按照转录方向进行反向互补处理，最终结果统一为启动子 5' 到 3' 方向。

---

### 6. Fielder promoter extraction / Fielder 启动子序列提取

Input Fielder gene IDs and retrieve Fielder promoter sequences.

输入 Fielder 基因号后，可批量获取 Fielder 启动子序列。

Fielder gene IDs usually follow this format:

Fielder 基因号通常为如下格式：

```text
TraesFLD5B01G105200
```

Current promoter definition:

当前启动子定义：

```text
2000 bp upstream of ATG
```

Features:

功能特点：

- Supports `TraesFLD...` gene IDs / 支持 `TraesFLD...` 格式基因号
- Fixed promoter length of 2000 bp / 固定提取 ATG 上游 2000 bp
- Strand-aware extraction / 考虑正负链方向
- Reverse-complement handling for negative-strand genes / 对负链基因进行反向互补处理
- FASTA output / 支持 FASTA 格式下载

---

### 7. JASPAR Plants promoter motif analysis / JASPAR Plants 启动子 motif 分析

Paste promoter FASTA records or plain DNA sequences and scan them against the local JASPAR CORE Plants non-redundant PWM database.

粘贴启动子 FASTA 或纯 DNA 序列后，可以使用本地 JASPAR CORE Plants non-redundant PWM 数据库扫描潜在转录因子结合位点。

Local motif files are stored in:

本地 motif 数据文件位于：

```text
data/motif_db/jaspar_plants/
├── JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt
├── JASPAR2026_CORE_plants_non-redundant_pfms_meme.txt
├── ultimate_metadata_table_CORE.tsv
├── jaspar_plants_pwm.json
├── background_cs_promoter_2000.json
├── background_fielder_promoter_2000.json
├── jaspar_background_thresholds_uniform.json
├── jaspar_background_thresholds_cs_promoter.json
└── jaspar_background_thresholds_fielder_promoter.json
```

The Streamlit page reads `jaspar_plants_pwm.json` and lightweight precomputed score threshold tables during analysis. It does not load full background score distributions by default.

Streamlit 页面运行分析时读取 `jaspar_plants_pwm.json` 和轻量级预计算 score threshold 表；默认不会加载完整背景分布 JSON。

Features:

功能特点：

- Supports FASTA and plain DNA input / 支持 FASTA 和纯 DNA 输入
- Supports multiple FASTA records / 支持多条 FASTA 序列
- Scans both forward and reverse-complement strands / 支持正链和反向互补链扫描
- Adjustable relative score pre-filter cutoff, default `0.90` / 可调 relative score 初筛阈值，默认 `0.90`
- Background score threshold p-level grading / 基于背景 score cutoff 的 p-level 分级
- Lightweight threshold JSON for Streamlit / 使用适合 Streamlit 的小型阈值 JSON
- Background options: uniform, Chinese Spring promoter, or Fielder promoter thresholds / 背景阈值表可选择均匀背景、中国春启动子背景或 Fielder 启动子背景
- Main table reports candidates passing the selected p-level threshold / 主结果表展示通过所选 p-level 阈值的候选结果
- Optional motif keyword filter by `matrix_id` or TF name / 可按 `matrix_id` 或 TF name 关键词筛选 motif
- Outputs detailed TF binding site candidates / 输出潜在 TF binding site 明细表
- Outputs motif-level summary table / 输出 motif 命中汇总表
- Supports CSV downloads / 支持 CSV 下载
- Uses a maximum hit limit to avoid freezing the web page / 使用最大 hits 数限制避免页面卡顿

Common output fields include:

常见输出字段包括：

```text
sequence_id
matrix_id
tf_name
consensus
motif_length
start
end
strand
matched_seq
raw_score
relative_score
p_level
confidence_level
significant
distance_to_sequence_end
species
tax_group
family
class
collection
```

Important interpretation note:

重要解释提醒：

PWM hits are predicted sequence matches only. The module uses precomputed background score cutoffs to assign p-level and confidence levels, which helps reduce false positives compared with using relative score alone while keeping Streamlit fast. This threshold mode does not output exact p-values or q-values. Significant-looking candidates still do not prove real TF binding or real regulatory relationships. Please combine expression data, conservation, ATAC-seq, ChIP-seq, or experimental validation.

PWM 命中只是序列层面的预测。本模块使用预计算背景 score cutoff 给结果分配 p-level 和 confidence level，相比单纯使用 relative score 可以减少假阳性，同时保持 Streamlit 快速运行。该 threshold 模式不输出精确 p-value 或 q-value。即使候选结果显著，也不等同于真实结合或真实调控证据。建议结合表达数据、保守性、ATAC-seq、ChIP-seq 或实验验证进一步确认。

---

#### Precomputed background thresholds / 预计算背景阈值表

For fast Streamlit analysis, precompute lightweight PWM background score thresholds offline:

为了让 Streamlit 页面快速分析，建议先离线预计算轻量级 PWM 背景 score cutoff 阈值表：

```bash
python scripts/build_jaspar_background.py
```

This generates `jaspar_background_thresholds_uniform.json` by default, and also generates Chinese Spring / Fielder promoter threshold files if their empirical background JSON files already exist.

该命令默认生成 `jaspar_background_thresholds_uniform.json`；如果中国春或 Fielder 启动子经验背景 JSON 已存在，也会同时生成对应阈值表。

To compute empirical promoter backgrounds from the local promoter SQLite databases first, run:

如需先从本地启动子 SQLite 数据库统计经验背景，运行：

```bash
python scripts/build_promoter_background.py
```

Then run:

然后运行：

```bash
python scripts/build_jaspar_background.py
```

The Streamlit page reports whether a lightweight threshold table has been loaded. If the selected threshold JSON is missing, the page asks you to run `python scripts/build_jaspar_background.py`.

Streamlit 页面会提示是否已加载轻量级阈值表。如果所选 threshold JSON 不存在，页面会提示先运行 `python scripts/build_jaspar_background.py`。

---

### 8. GO enrichment analysis / GO 富集分析

Input a DEG list with one gene ID per line and perform GO enrichment analysis.

输入差异基因列表，一行一个基因号，即可进行 GO 富集分析。

The GO enrichment module outputs:

GO 富集分析模块输出内容包括：

- Full GO enrichment result table / GO 富集结果总表
- Significant GO term table / 显著 GO 条目表
- GO enrichment bar plot / GO 富集条形图
- Analysis summary table / 分析摘要表

The enrichment analysis is based on local wheat GO annotation files and uses hypergeometric testing with multiple-testing correction.

GO 富集分析基于本地小麦 GO 注释文件，使用超几何检验进行富集分析，并进行多重检验校正。

The GO enrichment plot uses qvalue as the significance indicator. By default, it shows the top 15 most significant terms by qvalue for each GO ontology category (BP, CC, and MF), and this number can be adjusted on the page. Smaller qvalue means stronger enrichment significance. The plot style is kept close to the KEGG enrichment plots for easier comparison.

GO 富集图使用 qvalue 表示显著性，qvalue 越小表示富集越显著。默认情况下，每个 GO 大类（BP、CC、MF）展示 qvalue 最显著的前 15 个 term，可在页面中调整；图形风格已尽量与 KEGG 富集图保持一致，便于比较。

Local GO annotation files are stored in:

本地 GO 注释文件位于：

```text
data/go_mapping/
├── TERM2GENE_protein_coding.tsv
├── TERM2NAME_protein_coding.tsv
├── wheat_go_metadata.tsv
└── wheat_protein_coding_genes.tsv
```

File descriptions:

文件说明：

```text
TERM2GENE_protein_coding.tsv
```

GO term to gene mapping table.

GO 条目与基因之间的映射表。

```text
TERM2NAME_protein_coding.tsv
```

GO term ID to GO term name mapping table.

GO 编号与 GO 名称之间的映射表。

```text
wheat_go_metadata.tsv
```

GO annotation metadata.

GO 注释元数据信息。

```text
wheat_protein_coding_genes.tsv
```

Background gene set used for enrichment analysis.

用于富集分析的背景基因集。

GO enrichment workflow:

GO 富集分析流程：

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

Notes:

注意事项：

- GO enrichment results depend on the completeness of the local GO annotation.
- Genes without GO annotation will not contribute to GO enrichment testing.
- Enrichment results indicate over-representation of GO terms, not direct biological causality.
- Results should be interpreted together with biological knowledge and experimental context.

- GO 富集结果依赖本地 GO 注释文件的完整性。
- 没有 GO 注释的基因不会进入 GO 富集检验。
- 富集结果只能说明某些 GO 条目在输入基因集中显著偏多，不能直接证明因果关系。
- 结果应结合生物学背景和实验设计进行解释。

---

### 9. KEGG enrichment analysis / KEGG 富集分析

Input a DEG list with one gene ID per line and perform KEGG pathway enrichment analysis.

输入差异基因列表，一行一个基因号，即可进行 KEGG 通路富集分析。

The KEGG enrichment module outputs:

KEGG 富集分析模块输出内容包括：

- Full KEGG enrichment result table / KEGG 富集结果总表
- Significant pathway table / 显著通路结果表
- KEGG bubble plot / KEGG 富集气泡图
- KEGG bar plot / KEGG 富集条形图
- Analysis summary table / KEGG 分析摘要表

This module is based on local gene-KO and KO-pathway mapping files. It does not depend on online KEGG queries during analysis.

该模块基于本地 gene-KO 和 KO-pathway 映射文件进行分析，运行时不依赖在线 KEGG 查询。

Local KEGG files are stored in:

本地 KEGG 注释文件位于：

```text
data/kegg_mapping/
├── gene2ko_clean.tsv
├── kegg_ko2pathway.tsv
└── kegg_pathway2name.tsv
```

File descriptions:

文件说明：

```text
gene2ko_clean.tsv
```

Mapping between wheat genes and KEGG Orthology identifiers.

小麦基因与 KEGG Orthology，也就是 KO 编号之间的映射表。

Example:

示例：

```text
TraesCS1A02G000100    K00001
TraesCS1A02G000200    K14488
```

```text
kegg_ko2pathway.tsv
```

Mapping between KO identifiers and KEGG pathway IDs.

KO 编号与 KEGG pathway 编号之间的映射表。

Example:

示例：

```text
K00001    map00010
K14488    map04075
```

```text
kegg_pathway2name.tsv
```

Mapping between KEGG pathway IDs and pathway names.

KEGG pathway 编号与通路名称之间的映射表。

Example:

示例：

```text
map00010    Glycolysis / Gluconeogenesis
map04075    Plant hormone signal transduction
```

KEGG enrichment workflow:

KEGG 富集分析流程：

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

Common output fields:

常见输出字段说明：

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

Notes:

注意事项：

- KEGG enrichment is performed at the KO level rather than directly at the raw gene ID level.
- Genes without KO annotation will be excluded from KEGG enrichment testing.
- One gene may correspond to one or more KO identifiers.
- One KO identifier may participate in multiple KEGG pathways.
- KEGG enrichment indicates over-representation of pathways among input genes.
- KEGG enrichment does not directly indicate whether a pathway is activated or repressed.
- To infer pathway activation or repression, RNA-seq fold-change direction, expression pattern, and biological context should be considered together.

- KEGG 富集是在 KO 层面进行统计，而不是直接按原始基因号统计。
- 没有 KO 注释的基因不会进入 KEGG 富集检验。
- 一个基因可能对应一个或多个 KO 编号。
- 一个 KO 编号也可能参与多个 KEGG pathway。
- KEGG 富集结果说明某些通路在输入基因对应的 KO 中显著偏多。
- KEGG 富集结果不能直接说明通路被激活或被抑制。
- 如果需要判断通路上调或下调，需要结合 RNA-seq 的 log2FoldChange、表达趋势和具体生物学背景进一步解释。

---

## Database Design / 数据库设计

The original large SQLite database has been partitioned into multiple smaller query-ready SQLite databases.

原始大型 SQLite 数据库已经被拆分为多个可以直接查询的小型 SQLite 数据库。

This design avoids loading, merging, or decompressing a large database during web app startup.

这种设计避免了网站启动时加载、合并或解压大型数据库，从而更适合 GitHub 和 Streamlit Cloud 部署。

Current database structure:

当前数据库结构：

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

Example:

示例：

```text
TraesCS5B02G233300 → 5B
TraesFLD5B01G105200 → 5B
```

Then the program opens only the corresponding small SQLite database.

然后程序只打开对应染色体的小型 SQLite 数据库进行查询。

---

## Project Structure / 项目结构

```text
.
├── app.py
├── app_shared.py
├── requirements.txt
├── runtime.txt
├── start.bat
├── test.py
├── test_jaspar_scan.py
│
├── data/
│   ├── TipCode.jpg
│   ├── db/
│   ├── go_mapping/
│   ├── kegg_mapping/
│   └── motif_db/
│       └── jaspar_plants/
│           ├── JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt
│           ├── JASPAR2026_CORE_plants_non-redundant_pfms_meme.txt
│           ├── ultimate_metadata_table_CORE.tsv
│           ├── jaspar_plants_pwm.json
│           ├── background_cs_promoter_2000.json
│           ├── background_fielder_promoter_2000.json
│           ├── jaspar_background_thresholds_uniform.json
│           ├── jaspar_background_thresholds_cs_promoter.json
│           └── jaspar_background_thresholds_fielder_promoter.json
│
├── scripts/
│   ├── build_jaspar_background.py
│   ├── build_jaspar_pwm.py
│   ├── build_promoter_background.py
│   ├── split_sqlite_db.py
│   └── split_gene_structure_feature.py
│
├── sections/
│   ├── readme_page.py
│   ├── gene_info_page.py
│   ├── sequences_page.py
│   ├── homolog_page.py
│   ├── fielder_promoter_page.py
│   ├── cs_promoter_page.py
│   ├── motif_page.py
│   ├── go_page.py
│   └── kegg_page.py
│
└── utils/
    ├── db_query.py
    ├── go_enrichment.py
    ├── jaspar_pwm_scan.py
    └── kegg_enrichment.py
```

---

## Installation / 本地安装

Clone the repository:

克隆仓库：

```bash
git clone https://github.com/songsofdawn/wheat_genes_operation.git
cd wheat_genes_operation
```

Create a Python environment:

创建 Python 环境：

```bash
conda create -n wheattoolkit python=3.11 -y
conda activate wheattoolkit
```

Install dependencies:

安装依赖：

```bash
pip install -r requirements.txt
```

---

## Run Locally / 本地运行

Start the Streamlit app:

启动 Streamlit 应用：

```bash
streamlit run app.py
```

Then open the URL shown in the terminal.

然后打开终端中显示的网址。

Usually:

通常为：

```text
http://localhost:8501
```

---

## JASPAR PWM Scanner Test / JASPAR PWM 扫描测试

After `data/motif_db/jaspar_plants/jaspar_plants_pwm.json` is available, run:

当 `data/motif_db/jaspar_plants/jaspar_plants_pwm.json` 已经生成后，可以运行：

```bash
python test_jaspar_scan.py
```

The script checks whether the local JASPAR PWM JSON can be loaded, reads `jaspar_background_thresholds_uniform.json` if it exists, parses one test promoter, runs the threshold scanner, and writes:

该脚本会检查本地 JASPAR PWM JSON 能否读取；如果 `jaspar_background_thresholds_uniform.json` 存在，会优先读取它；然后解析一条测试启动子序列，调用 threshold 扫描函数，并输出：

```text
test_jaspar_scan_results.csv
```

This CSV is a local test artifact and is ignored by Git.

该 CSV 是本地测试产物，已加入 Git 忽略规则。

---

## Example Input / 示例输入

Chinese Spring gene IDs:

中国春基因号：

```text
TraesCS2D02G571200
TraesCS5B02G233300
TraesCS1A02G000100
```

Fielder gene IDs:

Fielder 基因号：

```text
TraesFLD5B01G105200
```

For batch analysis, use a TXT file with one gene ID per line.

批量分析时，请使用 TXT 文件，并保证每行一个基因号。

---

## Deployment on Streamlit Cloud / Streamlit Cloud 部署

This project can be deployed on Streamlit Cloud.

本项目可以部署到 Streamlit Cloud。

Deployment note: `data/db/` is required for the online app and is relatively large. On Streamlit Cloud, pay attention to repository size, cold-start time, disk usage, and memory limits. Do not commit local temporary files, cache folders, large intermediate JSON files, or local absolute paths.

部署提示：`data/db/` 是在线查询必需的数据目录，但体积较大。部署到 Streamlit Cloud 时需要注意仓库体积、冷启动时间、磁盘占用和内存限制。不要提交本地临时文件、缓存目录、大型中间 JSON 或本地绝对路径。

Recommended settings:

推荐设置：

```text
Repository: songsofdawn/wheat_genes_operation
Branch: main
Main file path: app.py
Python version: 3.11
```

The following files and directories should be included in the GitHub repository:

GitHub 仓库中应包含以下文件和目录：

```text
app.py
requirements.txt
runtime.txt
sections/
utils/
data/db/
data/go_mapping/
data/kegg_mapping/
data/motif_db/jaspar_plants/jaspar_plants_pwm.json
data/motif_db/jaspar_plants/background_cs_promoter_2000.json
data/motif_db/jaspar_plants/background_fielder_promoter_2000.json
data/motif_db/jaspar_plants/jaspar_background_thresholds_uniform.json
data/motif_db/jaspar_plants/jaspar_background_thresholds_cs_promoter.json
data/motif_db/jaspar_plants/jaspar_background_thresholds_fielder_promoter.json
data/TipCode.jpg
```

Do not upload the original large database:

不要上传原始大型数据库：

```text
wheat_toolkit.db
```

The online version uses the partitioned SQLite databases in `data/db/`.

在线版本使用 `data/db/` 中的分库 SQLite 数据库。

The JASPAR motif page requires `data/motif_db/jaspar_plants/jaspar_plants_pwm.json`. For fast Streamlit analysis, include the lightweight `jaspar_background_thresholds_*.json` files.

JASPAR motif 分析页面需要 `data/motif_db/jaspar_plants/jaspar_plants_pwm.json`。为了让 Streamlit 快速分析，建议同时包含轻量级 `jaspar_background_thresholds_*.json` 文件。

---

## Git Ignore Rules / Git 忽略规则

Recommended `.gitignore`:

推荐 `.gitignore`：

```gitignore
wheat_toolkit.db

*.db-journal
*.db-wal
*.db-shm

__pycache__/
*.pyc
*.pyo

.venv/
venv/
env/

.DS_Store
Thumbs.db

.streamlit/secrets.toml

test_jaspar_scan_results.csv
data/motif_db/jaspar_plants/jaspar_background_*_100000.json
data/motif_db/jaspar_plants/jaspar_background_*_50000.json
data/motif_db/jaspar_plants/jaspar_background_*_10000.json
data/motif_db/jaspar_plants/jaspar_background_*_cdf_*.json
```

Do not ignore:

不要忽略：

```text
data/db/
data/motif_db/jaspar_plants/jaspar_plants_pwm.json
```

because the small partitioned SQLite databases and the generated JASPAR PWM JSON are required for online deployment.

因为 `data/db/` 中的小型 SQLite 分库和生成好的 JASPAR PWM JSON 是在线部署所必需的。

---

## Database Construction Scripts / 数据库构建脚本

The `scripts/` directory contains scripts used to split the original SQLite database into smaller query-ready databases.

`scripts/` 目录包含用于拆分原始 SQLite 数据库的脚本。

```text
scripts/
├── build_jaspar_background.py
├── build_jaspar_pwm.py
├── build_promoter_background.py
├── split_sqlite_db.py
└── split_gene_structure_feature.py
```

These scripts are mainly used for database and motif resource construction and maintenance. `build_promoter_background.py` computes empirical A/C/G/T backgrounds from local promoter SQLite databases. `build_jaspar_background.py` precomputes lightweight PWM score thresholds for Streamlit.

这些脚本主要用于数据库和 motif 资源构建维护。`build_promoter_background.py` 会从本地启动子 SQLite 数据库统计 A/C/G/T 经验背景；`build_jaspar_background.py` 会预计算轻量级 PWM score cutoff 阈值表，供 Streamlit 页面快速分级。

---

## Notes / 注意事项

- This tool is designed for wheat gene list processing and functional genomics analysis.
- The current promoter definition is fixed as 2000 bp upstream of ATG.
- Homolog relationships are computationally inferred and should be interpreted with caution.
- JASPAR PWM motif hits are potential TF binding site predictions, not direct evidence of real binding or regulation.
- GO enrichment results depend on the quality and completeness of local GO annotation files.
- KEGG enrichment results depend on the quality and completeness of local gene-KO and KO-pathway mapping files.
- Enrichment results should be interpreted together with biological knowledge and experimental validation.

- 本工具主要用于小麦基因列表处理和功能基因组学分析。
- 当前启动子定义固定为 ATG 上游 2000 bp。
- 同源关系为计算推断结果，应谨慎解释。
- JASPAR PWM motif 命中是潜在 TF binding site 预测，不等同于真实结合或真实调控证据。
- GO 富集结果依赖本地 GO 注释文件的质量和完整性。
- KEGG 富集结果依赖本地 gene-KO 和 KO-pathway 映射文件的质量和完整性。
- 富集分析结果应结合生物学知识和实验验证进行综合判断。

---

## Future Development / 后续开发计划

Potential future modules include:

未来可扩展模块包括：

- Gene structure visualization / 基因结构可视化
- Expression profile visualization / 表达谱可视化
- Motif result visualization / motif 结果可视化
- Co-expression network query / 共表达网络查询
- REST API service / API 接口服务
- Batch result packaging and ZIP download / 批量结果打包下载

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

The developer is not responsible for incorrect biological conclusions caused by inappropriate use of the tool or unverified downstream interpretation.

对于因不当使用本工具或未经验证的下游解释所导致的错误生物学结论，开发者不承担责任。
