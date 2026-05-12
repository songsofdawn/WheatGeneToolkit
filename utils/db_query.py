# -*- coding: utf-8 -*-

import json
import re
import sqlite3
from functools import lru_cache
from pathlib import Path
from typing import Optional, Any

import pandas as pd

try:
    import streamlit as st
except Exception:  # pragma: no cover - 兼容非 Streamlit 脚本环境
    st = None


# ============================================================
# 项目路径与分库 manifest
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DB_DIR = PROJECT_ROOT / "data" / "db"
MANIFEST_FILE = DB_DIR / "manifest.json"
DATABASE_MISSING_MESSAGE = "数据库文件缺失，请检查 data/db/manifest.json 和分库 SQLite 文件是否存在。"


class DatabaseUnavailableError(RuntimeError):
    """数据库资源不可用时抛出的友好异常，避免模块导入阶段直接崩溃。"""


def _get_file_cache_key(path: Path) -> tuple[str, int, int]:
    """
    Build a cache key that invalidates automatically when a file changes.
    """
    resolved = path.resolve()
    stat = resolved.stat()
    return str(resolved), stat.st_mtime_ns, stat.st_size


@lru_cache(maxsize=1)
def _load_manifest_cached(path_str: str, _mtime_ns: int, _size: int) -> dict:
    with Path(path_str).open("r", encoding="utf-8") as f:
        return json.load(f)


def _load_manifest() -> dict:
    """
    读取 data/db/manifest.json。
    """
    if not MANIFEST_FILE.exists():
        raise FileNotFoundError(
            f"找不到分库清单文件: {MANIFEST_FILE}\n"
            f"请确认你已经完成数据库拆分，并且 data/db/manifest.json 存在。"
        )

    return _load_manifest_cached(*_get_file_cache_key(MANIFEST_FILE))


MANIFEST = None


def _get_manifest() -> dict:
    """
    懒加载 data/db/manifest.json。
    这样即使云端部署时数据库文件缺失，Streamlit 入口和 ReadMe 页面也能正常打开。
    """
    global MANIFEST

    if MANIFEST is not None:
        return MANIFEST

    try:
        MANIFEST = _load_manifest()
        return MANIFEST
    except FileNotFoundError as exc:
        raise DatabaseUnavailableError(f"{DATABASE_MISSING_MESSAGE}\n{exc}") from exc
    except Exception as exc:
        raise DatabaseUnavailableError(f"{DATABASE_MISSING_MESSAGE}\n读取 manifest 失败: {exc}") from exc


# ============================================================
# 基础工具函数
# ============================================================

def _quote_ident(name: str) -> str:
    """
    SQLite 字段/表名安全引用。
    """
    return '"' + str(name).replace('"', '""') + '"'


def _normalize_input_id(input_id: Any) -> str:
    """
    清理用户输入的基因号。
    """
    if input_id is None:
        return ""

    s = str(input_id).strip()

    # 如果用户从 FASTA header 复制了 >gene_id|xxx，只取第一个字段
    s = s.lstrip(">").split("|")[0].strip()

    return s


def _remove_transcript_suffix(gene_id: str) -> str:
    """
    把 TraesCSxxx.1 这类 transcript 后缀去掉。
    例如：
    TraesCS5B02G233300.1 -> TraesCS5B02G233300
    """
    gene_id = _normalize_input_id(gene_id)
    return re.sub(r"\.\d+$", "", gene_id)


def parse_chr_from_gene_id(gene_id: Any) -> str:
    """
    根据基因号解析染色体编号。

    支持：
    TraesCS5B02G233300
    TraesCS5B02G233300.1
    TraesFLD5B01G105200
    """
    s = _normalize_input_id(gene_id)

    m = re.search(r"TraesCS([1-7][ABD])", s)
    if m:
        return m.group(1)

    m = re.search(r"TraesFLD([1-7][ABD])", s)
    if m:
        return m.group(1)

    return "unknown"


def _connect_db(db_path: Path) -> sqlite3.Connection:
    """
    获取 SQLite 连接。
    """
    if not db_path.exists():
        raise FileNotFoundError(f"找不到数据库文件: {db_path}")

    return get_sqlite_connection(*_get_file_cache_key(db_path))


def _open_readonly_connection(path_str: str, _mtime_ns: int, _size: int) -> sqlite3.Connection:
    """
    以只读 URI 打开 SQLite，并复用连接，减少 Streamlit Cloud 上的反复 IO。
    """
    db_path = Path(path_str)
    if not db_path.exists():
        raise FileNotFoundError(f"找不到数据库文件: {db_path}")

    uri = db_path.resolve().as_uri() + "?mode=ro"
    conn = sqlite3.connect(uri, uri=True, check_same_thread=False)
    conn.execute("PRAGMA query_only = ON")
    conn.execute("PRAGMA temp_store = MEMORY")
    conn.execute("PRAGMA cache_size = -20000")
    return conn


if st is not None:
    get_sqlite_connection = st.cache_resource(show_spinner=False)(_open_readonly_connection)
else:
    get_sqlite_connection = lru_cache(maxsize=128)(_open_readonly_connection)


@lru_cache(maxsize=2048)
def _read_sql_query_cached(
    db_path_str: str,
    _mtime_ns: int,
    _size: int,
    sql: str,
    params: tuple,
) -> pd.DataFrame:
    conn = _connect_db(Path(db_path_str))
    return pd.read_sql_query(sql, conn, params=params)


def _read_sql_query_uncached(db_path: Path, sql: str, params: tuple = ()) -> pd.DataFrame:
    conn = _connect_db(db_path)
    return pd.read_sql_query(sql, conn, params=params)


def get_connection(db_file: Optional[str] = None) -> sqlite3.Connection:
    """
    兼容旧代码的连接函数。

    注意：
    - 旧版本默认连接 wheat_toolkit.db；
    - 新版本不再推荐直接调用这个函数；
    - 如果 db_file 不传，则默认连接 gene_core.db，仅用于少数兼容场景。
    """
    if db_file:
        return _connect_db(PROJECT_ROOT / db_file)

    db_path = _get_table_db_path("gene_core")
    return _connect_db(db_path)


def _get_table_info(table: str) -> dict:
    """
    从 manifest 中获得某个表的拆分信息。
    """
    manifest = _get_manifest()
    try:
        return manifest["tables"][table]
    except KeyError:
        raise KeyError(
            f"manifest.json 中找不到表 {table}。\n"
            f"请检查 data/db/manifest.json 是否包含该表。"
        )


def _get_table_db_path(table: str, gene_id: Optional[str] = None) -> Path:
    """
    根据表名和 gene_id 返回对应的小 SQLite 数据库路径。

    whole_table:
        data/db/core/gene_core.db

    by_chr:
        data/db/gene_promoter_sequence/gene_promoter_sequence_5B.db
    """
    info = _get_table_info(table)
    strategy = info.get("strategy")

    if strategy == "whole_table":
        return PROJECT_ROOT / info["path"]

    if strategy in {"by_chr", "by_chr_from_gene_id"}:
        if gene_id is None:
            raise ValueError(f"表 {table} 是按染色体拆分的，必须提供 gene_id。")

        chr_name = parse_chr_from_gene_id(gene_id)
        shards = info.get("shards", {})

        # 兼容 shards 是 dict 的情况：
        # "shards": {"5B": {"path": "..."}}
        if isinstance(shards, dict):
            if chr_name in shards:
                return PROJECT_ROOT / shards[chr_name]["path"]

            if "unknown" in shards:
                return PROJECT_ROOT / shards["unknown"]["path"]

        # 兼容 shards 是 list 的情况：
        # "shards": [{"chr": "5B", "path": "..."}]
        if isinstance(shards, list):
            for shard in shards:
                if shard.get("chr") == chr_name:
                    return PROJECT_ROOT / shard["path"]

            for shard in shards:
                if shard.get("chr") == "unknown":
                    return PROJECT_ROOT / shard["path"]

        raise FileNotFoundError(
            f"找不到表 {table} 对应染色体 {chr_name} 的分库。\n"
            f"gene_id = {gene_id}"
        )

    raise ValueError(f"未知拆分策略: table={table}, strategy={strategy}")


def _fetch_df(
    table: str,
    sql: str,
    params: tuple = (),
    gene_id_for_route: Optional[str] = None,
) -> pd.DataFrame:
    """
    通用查询函数，自动根据 table/gene_id 路由到正确的小数据库。
    """
    try:
        db_path = _get_table_db_path(table, gene_id_for_route)
        return _read_sql_query_cached(
            *_get_file_cache_key(db_path),
            sql,
            tuple(params),
        ).copy()
    except FileNotFoundError as exc:
        raise DatabaseUnavailableError(f"{DATABASE_MISSING_MESSAGE}\n缺失文件: {exc}") from exc


def _chunked(values, chunk_size=800):
    values = list(values)
    for i in range(0, len(values), chunk_size):
        yield values[i : i + chunk_size]


def _placeholders(count: int) -> str:
    return ",".join(["?"] * count)


def _unique_preserve_order(values):
    seen = set()
    result = []
    for value in values:
        if value is None:
            continue
        value = str(value).strip()
        if not value or value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def _fetch_in_from_db_path(
    db_path: Path,
    table: str,
    column: str,
    values,
    order_by: Optional[str] = None,
) -> pd.DataFrame:
    values = _unique_preserve_order(values)
    if not values:
        return pd.DataFrame()

    frames = []
    quoted_table = _quote_ident(table)
    quoted_col = _quote_ident(column)
    for chunk in _chunked(values):
        sql = f"""
        SELECT *
        FROM {quoted_table}
        WHERE {quoted_col} IN ({_placeholders(len(chunk))})
        """
        if order_by:
            sql += f"\nORDER BY {order_by}"
        frames.append(_read_sql_query_uncached(db_path, sql, tuple(chunk)))

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _fetch_in_whole_table(
    table: str,
    column: str,
    values,
    order_by: Optional[str] = None,
) -> pd.DataFrame:
    db_path = _get_table_db_path(table)
    return _fetch_in_from_db_path(db_path, table, column, values, order_by=order_by)


def _fetch_in_sharded_table(
    table: str,
    column: str,
    values,
    order_by: Optional[str] = None,
) -> pd.DataFrame:
    groups = {}
    for value in _unique_preserve_order(values):
        db_path = _get_table_db_path(table, value)
        groups.setdefault(db_path, []).append(value)

    frames = [
        _fetch_in_from_db_path(db_path, table, column, group_values, order_by=order_by)
        for db_path, group_values in groups.items()
    ]
    frames = [df for df in frames if df is not None and not df.empty]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def _add_input_order(df: pd.DataFrame, gene_ids, id_map: dict, result_id_col="primary_gene_id") -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()

    order_rows = []
    for idx, input_gene_id in enumerate(gene_ids):
        resolved_id = id_map.get(input_gene_id)
        if resolved_id:
            order_rows.append({
                "input_gene_id": input_gene_id,
                "input_order": idx,
                result_id_col: resolved_id,
            })

    if not order_rows:
        return pd.DataFrame()

    order_df = pd.DataFrame(order_rows)
    merged = order_df.merge(df, on=result_id_col, how="inner")
    return merged.sort_values("input_order").reset_index(drop=True)


# ============================================================
# ID 转换
# ============================================================

def get_primary_gene_id(input_id: str, db_file: Optional[str] = None) -> Optional[str]:
    """
    根据输入的二代号 / 三代号 / alias，返回统一主键 primary_gene_id。

    db_file 参数保留是为了兼容旧 app.py，不再实际使用。
    """
    input_id = _normalize_input_id(input_id)
    if not input_id:
        return None

    input_id_no_suffix = _remove_transcript_suffix(input_id)

    candidates = []
    for x in [input_id, input_id_no_suffix]:
        if x and x not in candidates:
            candidates.append(x)

    # 1. 先查 alias 表
    for x in candidates:
        sql = """
        SELECT primary_gene_id
        FROM gene_alias
        WHERE alias_value = ?
        LIMIT 1
        """
        df = _fetch_df("gene_alias", sql, (x,))
        if not df.empty:
            return df.loc[0, "primary_gene_id"]

    # 2. 再查 gene_core 主键
    for x in candidates:
        sql = """
        SELECT primary_gene_id
        FROM gene_core
        WHERE primary_gene_id = ?
        LIMIT 1
        """
        df = _fetch_df("gene_core", sql, (x,))
        if not df.empty:
            return df.loc[0, "primary_gene_id"]

    # 3. 再查三代号 gene_id_v3
    for x in candidates:
        sql = """
        SELECT primary_gene_id
        FROM gene_core
        WHERE gene_id_v3 = ?
        LIMIT 1
        """
        df = _fetch_df("gene_core", sql, (x,))
        if not df.empty:
            return df.loc[0, "primary_gene_id"]

    return None


def _candidate_ids(input_id: Any) -> list[str]:
    input_id = _normalize_input_id(input_id)
    input_id_no_suffix = _remove_transcript_suffix(input_id)
    return _unique_preserve_order([input_id, input_id_no_suffix])


@lru_cache(maxsize=512)
def _resolve_primary_gene_ids_many_cached(gene_ids_tuple: tuple) -> tuple:
    gene_ids = [_normalize_input_id(gene_id) for gene_id in gene_ids_tuple if _normalize_input_id(gene_id)]
    gene_ids = _unique_preserve_order(gene_ids)
    candidate_by_input = {gene_id: _candidate_ids(gene_id) for gene_id in gene_ids}
    all_candidates = _unique_preserve_order(
        candidate for candidates in candidate_by_input.values() for candidate in candidates
    )

    alias_map = {}
    primary_map = {}
    v3_map = {}

    if all_candidates:
        alias_df = _fetch_in_whole_table("gene_alias", "alias_value", all_candidates)
        if not alias_df.empty:
            alias_map = dict(zip(alias_df["alias_value"], alias_df["primary_gene_id"]))

        core_primary_df = _fetch_in_whole_table("gene_core", "primary_gene_id", all_candidates)
        if not core_primary_df.empty:
            primary_map = dict(zip(core_primary_df["primary_gene_id"], core_primary_df["primary_gene_id"]))

        core_v3_df = _fetch_in_whole_table("gene_core", "gene_id_v3", all_candidates)
        if not core_v3_df.empty:
            v3_map = dict(zip(core_v3_df["gene_id_v3"], core_v3_df["primary_gene_id"]))

    rows = []
    for gene_id in gene_ids:
        resolved = None
        candidates = candidate_by_input.get(gene_id, [])
        for candidate in candidates:
            if candidate in alias_map:
                resolved = alias_map[candidate]
                break
        if resolved is None:
            for candidate in candidates:
                if candidate in primary_map:
                    resolved = primary_map[candidate]
                    break
        if resolved is None:
            for candidate in candidates:
                if candidate in v3_map:
                    resolved = v3_map[candidate]
                    break
        rows.append((gene_id, resolved))

    return tuple(rows)


def resolve_primary_gene_ids_many(gene_ids) -> dict:
    """
    批量解析输入基因号到 primary_gene_id，保留输入顺序并缓存映射结果。
    """
    cleaned = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    return dict(_resolve_primary_gene_ids_many_cached(tuple(cleaned)))


def get_gene_core_many(gene_ids) -> pd.DataFrame:
    gene_ids = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
    primary_ids = _unique_preserve_order(id_map.values())
    df = _fetch_in_whole_table("gene_core", "primary_gene_id", primary_ids)
    return _add_input_order(df, gene_ids, id_map, result_id_col="primary_gene_id")


def get_gene_annotations_many(gene_ids) -> pd.DataFrame:
    gene_ids = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
    primary_ids = _unique_preserve_order(id_map.values())
    df = _fetch_in_whole_table("gene_annotation", "primary_gene_id", primary_ids)
    return _add_input_order(df, gene_ids, id_map, result_id_col="primary_gene_id")


def get_transcript_core_many(gene_ids) -> pd.DataFrame:
    gene_ids = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
    primary_ids = _unique_preserve_order(id_map.values())
    df = _fetch_in_whole_table(
        "transcript_core",
        "primary_gene_id",
        primary_ids,
        order_by="primary_gene_id, is_canonical DESC, transcript_id ASC",
    )
    return _add_input_order(df, gene_ids, id_map, result_id_col="primary_gene_id")


def get_gene_sequence_resources_many(gene_ids) -> pd.DataFrame:
    gene_ids = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
    primary_ids = _unique_preserve_order(id_map.values())
    df = _fetch_in_sharded_table(
        "gene_sequence_resource",
        "primary_gene_id",
        primary_ids,
        order_by="primary_gene_id, transcript_id, sequence_type",
    )
    return _add_input_order(df, gene_ids, id_map, result_id_col="primary_gene_id")


def get_gene_promoters_many(gene_ids, genome="CS") -> pd.DataFrame:
    gene_ids = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    if str(genome).upper() == "FIELDER":
        return get_fielder_promoters_many(gene_ids)

    id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
    primary_ids = _unique_preserve_order(id_map.values())
    df = _fetch_in_sharded_table("gene_promoter_sequence", "primary_gene_id", primary_ids)
    return _add_input_order(df, gene_ids, id_map, result_id_col="primary_gene_id")


def get_fielder_promoters_many(gene_ids) -> pd.DataFrame:
    gene_ids = _unique_preserve_order(_remove_transcript_suffix(gene_id) for gene_id in gene_ids)
    id_map = {gene_id: gene_id for gene_id in gene_ids if gene_id}
    df = _fetch_in_sharded_table("fielder_promoter_sequence", "primary_gene_id", gene_ids)
    return _add_input_order(df, gene_ids, id_map, result_id_col="primary_gene_id")


def _add_homolog_input_order(df: pd.DataFrame, gene_ids, id_map: dict) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()
    order_df = pd.DataFrame(
        [
            {"input_gene_id": input_gene_id, "input_order": idx, "cs_gene_id": primary_id}
            for idx, input_gene_id in enumerate(gene_ids)
            for primary_id in [id_map.get(input_gene_id)]
            if primary_id
        ]
    )
    if order_df.empty:
        return pd.DataFrame()
    return order_df.merge(df, on="cs_gene_id", how="inner").sort_values(
        ["input_order", "rank_within_cs"]
    ).reset_index(drop=True)


def get_homologs_many(gene_ids) -> dict:
    gene_ids = _unique_preserve_order(_normalize_input_id(gene_id) for gene_id in gene_ids)
    id_map = resolve_primary_gene_ids_many(tuple(gene_ids))
    primary_ids = _unique_preserve_order(id_map.values())

    self_all = _fetch_in_whole_table(
        "cs_self_homolog_map",
        "cs_gene_id",
        primary_ids,
        order_by="cs_gene_id, CAST(rank_within_cs AS INTEGER), CAST(priority_score AS REAL) DESC",
    )
    fielder_all = _fetch_in_whole_table(
        "homolog_map",
        "cs_gene_id",
        primary_ids,
        order_by="cs_gene_id, CAST(rank_within_cs AS INTEGER), CAST(priority_score AS REAL) DESC",
    )

    self_all = _add_homolog_input_order(self_all, gene_ids, id_map)
    fielder_all = _add_homolog_input_order(fielder_all, gene_ids, id_map)

    self_best = (
        self_all[self_all["is_best_hit"].astype(str) == "1"].copy()
        if not self_all.empty and "is_best_hit" in self_all.columns
        else pd.DataFrame()
    )
    fielder_best = (
        fielder_all[fielder_all["is_best_hit"].astype(str) == "1"].copy()
        if not fielder_all.empty and "is_best_hit" in fielder_all.columns
        else pd.DataFrame()
    )

    found_primary_ids = set()
    for df in [self_best, fielder_best]:
        if not df.empty and "cs_gene_id" in df.columns:
            found_primary_ids.update(df["cs_gene_id"].dropna().astype(str))

    missing = [
        gene_id
        for gene_id in gene_ids
        if not id_map.get(gene_id) or id_map.get(gene_id) not in found_primary_ids
    ]

    return {
        "id_map": id_map,
        "self_best": self_best,
        "self_all": self_all,
        "fielder_best": fielder_best,
        "fielder_all": fielder_all,
        "missing": missing,
    }


# ============================================================
# 中国春 gene core / annotation / transcript
# ============================================================

def get_gene_core(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取 gene_core 基础信息。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM gene_core
    WHERE primary_gene_id = ?
    """
    return _fetch_df("gene_core", sql, (primary_gene_id,))


def get_gene_annotation(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取基因功能注释。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM gene_annotation
    WHERE primary_gene_id = ?
    """
    return _fetch_df("gene_annotation", sql, (primary_gene_id,))


def get_transcripts(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取某个基因对应的 transcript 列表。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM transcript_core
    WHERE primary_gene_id = ?
    ORDER BY is_canonical DESC, transcript_id ASC
    """
    return _fetch_df("transcript_core", sql, (primary_gene_id,))


# ============================================================
# 序列 / 启动子 / 结构
# ============================================================

def get_sequences(
    input_id: str,
    sequence_type: Optional[str] = None,
    db_file: Optional[str] = None,
) -> pd.DataFrame:
    """
    获取序列资源。

    sequence_type 可选:
        - cdna
        - cds
        - protein

    不传则返回全部。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    if sequence_type is None:
        sql = """
        SELECT *
        FROM gene_sequence_resource
        WHERE primary_gene_id = ?
        ORDER BY transcript_id, sequence_type
        """
        return _fetch_df(
            "gene_sequence_resource",
            sql,
            (primary_gene_id,),
            gene_id_for_route=primary_gene_id,
        )

    sql = """
    SELECT *
    FROM gene_sequence_resource
    WHERE primary_gene_id = ?
      AND sequence_type = ?
    ORDER BY transcript_id
    """
    return _fetch_df(
        "gene_sequence_resource",
        sql,
        (primary_gene_id, sequence_type),
        gene_id_for_route=primary_gene_id,
    )


def get_promoter(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取中国春启动子序列。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM gene_promoter_sequence
    WHERE primary_gene_id = ?
    """
    return _fetch_df(
        "gene_promoter_sequence",
        sql,
        (primary_gene_id,),
        gene_id_for_route=primary_gene_id,
    )


def get_gene_structure(
    input_id: str,
    transcript_id: Optional[str] = None,
    db_file: Optional[str] = None,
) -> pd.DataFrame:
    """
    获取基因结构信息。

    transcript_id 不传则返回该基因全部 transcript 的结构。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    if transcript_id is None:
        sql = """
        SELECT *
        FROM gene_structure_feature
        WHERE primary_gene_id = ?
        ORDER BY transcript_id, feature_order, start
        """
        return _fetch_df(
            "gene_structure_feature",
            sql,
            (primary_gene_id,),
            gene_id_for_route=primary_gene_id,
        )

    sql = """
    SELECT *
    FROM gene_structure_feature
    WHERE primary_gene_id = ?
      AND transcript_id = ?
    ORDER BY feature_order, start
    """
    return _fetch_df(
        "gene_structure_feature",
        sql,
        (primary_gene_id, transcript_id),
        gene_id_for_route=primary_gene_id,
    )


# ============================================================
# alias 搜索
# ============================================================

def search_alias(
    keyword: str,
    limit: int = 20,
    db_file: Optional[str] = None,
) -> pd.DataFrame:
    """
    模糊搜索基因号。
    """
    keyword = _normalize_input_id(keyword)
    if not keyword:
        return pd.DataFrame()

    keyword_like = f"%{keyword}%"
    limit = max(1, min(int(limit), 500))

    sql = f"""
    SELECT *
    FROM gene_alias
    WHERE alias_value LIKE ?
    LIMIT {limit}
    """
    return _fetch_df("gene_alias", sql, (keyword_like,))


# ============================================================
# bundle
# ============================================================

def get_full_gene_bundle(input_id: str, db_file: Optional[str] = None) -> dict:
    """
    一次性返回一个基因的常用信息。
    """
    primary_gene_id = get_primary_gene_id(input_id)

    if primary_gene_id is None:
        return {
            "primary_gene_id": None,
            "gene_core": pd.DataFrame(),
            "gene_annotation": pd.DataFrame(),
            "transcripts": pd.DataFrame(),
            "sequences": pd.DataFrame(),
            "promoter": pd.DataFrame(),
            "structure": pd.DataFrame(),
        }

    return {
        "primary_gene_id": primary_gene_id,
        "gene_core": get_gene_core(primary_gene_id),
        "gene_annotation": get_gene_annotation(primary_gene_id),
        "transcripts": get_transcripts(primary_gene_id),
        "sequences": get_sequences(primary_gene_id, None),
        "promoter": get_promoter(primary_gene_id),
        "structure": get_gene_structure(primary_gene_id, None),
    }


# ============================================================
# 中国春 -> Fielder 同源
# ============================================================

def get_fielder_best_hit(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取某个中国春基因最可信的 Fielder 对应基因。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM homolog_map
    WHERE cs_gene_id = ?
      AND is_best_hit = '1'
    ORDER BY CAST(rank_within_cs AS INTEGER)
    """
    return _fetch_df("homolog_map", sql, (primary_gene_id,))


def get_fielder_all_hits(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取某个中国春基因所有 Fielder 候选同源基因。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM homolog_map
    WHERE cs_gene_id = ?
    ORDER BY CAST(rank_within_cs AS INTEGER), CAST(priority_score AS REAL) DESC
    """
    return _fetch_df("homolog_map", sql, (primary_gene_id,))


# ============================================================
# Fielder 基因信息 / 启动子
# ============================================================

def get_fielder_gene_core(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取 Fielder gene_core 基础信息。

    input_id 直接是 Fielder 基因号，例如：
    TraesFLD1A01G000100
    """
    input_id = _remove_transcript_suffix(input_id)

    sql = """
    SELECT *
    FROM fielder_gene_core
    WHERE primary_gene_id = ?
    """
    return _fetch_df("fielder_gene_core", sql, (input_id,))


def get_fielder_promoter(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取 Fielder 启动子序列。

    input_id 直接是 Fielder 基因号，例如：
    TraesFLD1A01G000100
    """
    input_id = _remove_transcript_suffix(input_id)

    if not input_id:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM fielder_promoter_sequence
    WHERE primary_gene_id = ?
    """
    return _fetch_df(
        "fielder_promoter_sequence",
        sql,
        (input_id,),
        gene_id_for_route=input_id,
    )


# ============================================================
# 中国春自身同源
# ============================================================

def get_cs_self_best_hit(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取某个中国春基因最可信的自身同源基因。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM cs_self_homolog_map
    WHERE cs_gene_id = ?
      AND is_best_hit = '1'
    ORDER BY CAST(rank_within_cs AS INTEGER)
    """
    return _fetch_df("cs_self_homolog_map", sql, (primary_gene_id,))


def get_cs_self_all_hits(input_id: str, db_file: Optional[str] = None) -> pd.DataFrame:
    """
    获取某个中国春基因所有自身同源候选基因。
    """
    primary_gene_id = get_primary_gene_id(input_id)
    if primary_gene_id is None:
        return pd.DataFrame()

    sql = """
    SELECT *
    FROM cs_self_homolog_map
    WHERE cs_gene_id = ?
    ORDER BY CAST(rank_within_cs AS INTEGER), CAST(priority_score AS REAL) DESC
    """
    return _fetch_df("cs_self_homolog_map", sql, (primary_gene_id,))


# ============================================================
# 部署前自检
# ============================================================

def check_database_status() -> pd.DataFrame:
    """
    检查 manifest 中记录的所有数据库文件是否存在。
    """
    rows = []

    manifest = _get_manifest()

    for table, info in manifest.get("tables", {}).items():
        strategy = info.get("strategy")

        if strategy == "whole_table":
            path = PROJECT_ROOT / info["path"]
            rows.append({
                "table": table,
                "strategy": strategy,
                "shard": "",
                "path": str(path),
                "exists": path.exists(),
                "size_mb": round(path.stat().st_size / 1024 / 1024, 2) if path.exists() else None,
            })

        elif strategy in {"by_chr", "by_chr_from_gene_id"}:
            shards = info.get("shards", {})

            if isinstance(shards, dict):
                iterable = shards.items()
                for shard_name, shard_info in iterable:
                    path = PROJECT_ROOT / shard_info["path"]
                    rows.append({
                        "table": table,
                        "strategy": strategy,
                        "shard": shard_name,
                        "path": str(path),
                        "exists": path.exists(),
                        "size_mb": round(path.stat().st_size / 1024 / 1024, 2) if path.exists() else None,
                    })

            elif isinstance(shards, list):
                for shard_info in shards:
                    shard_name = shard_info.get("chr", "")
                    path = PROJECT_ROOT / shard_info["path"]
                    rows.append({
                        "table": table,
                        "strategy": strategy,
                        "shard": shard_name,
                        "path": str(path),
                        "exists": path.exists(),
                        "size_mb": round(path.stat().st_size / 1024 / 1024, 2) if path.exists() else None,
                    })

    return pd.DataFrame(rows)


if __name__ == "__main__":
    test_id = "TraesCS2D02G571200"

    print("=== 数据库文件检查 ===")
    status_df = check_database_status()
    print(status_df.head())
    print("缺失文件数:", (~status_df["exists"]).sum() if not status_df.empty else "未知")

    print("\n=== primary_gene_id ===")
    print(get_primary_gene_id(test_id))

    print("\n=== gene_core ===")
    print(get_gene_core(test_id).head())

    print("\n=== gene_annotation ===")
    print(get_gene_annotation(test_id).head())

    print("\n=== transcripts ===")
    print(get_transcripts(test_id).head())

    print("\n=== sequences ===")
    print(get_sequences(test_id).head())

    print("\n=== promoter ===")
    print(get_promoter(test_id).head())

    print("\n=== structure ===")
    print(get_gene_structure(test_id).head())

    print("\n=== fielder best hit ===")
    print(get_fielder_best_hit(test_id).head())

    print("\n=== cs self best hit ===")
    print(get_cs_self_best_hit(test_id).head())
