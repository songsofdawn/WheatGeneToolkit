# -*- coding: utf-8 -*-

"""
将 wheat_toolkit.db 拆分成多个可直接查询的小 SQLite 数据库。

目标结构：

data/db/
├── gene_sequence_resource/
│   ├── gene_sequence_resource_1A.db
│   ├── gene_sequence_resource_1B.db
│   └── ...
├── gene_promoter_sequence/
│   ├── gene_promoter_sequence_1A.db
│   ├── gene_promoter_sequence_1B.db
│   └── ...
├── fielder_promoter_sequence/
│   ├── fielder_promoter_sequence_1A.db
│   ├── fielder_promoter_sequence_1B.db
│   └── ...
├── core/
├── homolog/
└── manifest.json

运行方式：

cd <project-root>
python scripts\\split_sqlite_db.py
"""

import json
import re
import shutil
import sqlite3
import time
from pathlib import Path


PROJECT_ROOT = Path.cwd()
SOURCE_DB = PROJECT_ROOT / "wheat_toolkit.db"
OUT_DIR = PROJECT_ROOT / "data" / "db"

MAX_FILE_MB = 95
BATCH_SIZE = 5000


# 这些表按照染色体拆分
SPLIT_BY_CHR_TABLES = {
    "gene_sequence_resource",
    "gene_promoter_sequence",
    "fielder_promoter_sequence",
}

# 这些表作为核心表，通常不大，按表单独存为一个 db
CORE_TABLES = {
    "gene_core",
    "fielder_gene_core",
    "gene_alias",
    "gene_annotation",
    "gene_structure_feature",
    "transcript_core",
}

# 同源表，默认先整体复制；如果超过 95 MB，脚本会提示你后续再拆
HOMOLOG_TABLES = {
    "homolog_map",
    "cs_self_homolog_map",
}


def log(msg: str):
    print(msg, flush=True)


def mb(size_bytes: int) -> float:
    return size_bytes / 1024 / 1024


def file_size_mb(path: Path) -> float:
    if not path.exists():
        return 0.0
    return mb(path.stat().st_size)


def relpath(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")


def connect_sqlite(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    return conn


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def get_tables(conn: sqlite3.Connection):
    rows = conn.execute("""
        SELECT name
        FROM sqlite_master
        WHERE type='table'
          AND name NOT LIKE 'sqlite_%'
        ORDER BY name
    """).fetchall()
    return [r[0] for r in rows]


def get_columns(conn: sqlite3.Connection, table: str):
    rows = conn.execute(f"PRAGMA table_info({quote_ident(table)})").fetchall()
    return [r[1] for r in rows]


def get_create_table_sql(conn: sqlite3.Connection, table: str) -> str:
    row = conn.execute("""
        SELECT sql
        FROM sqlite_master
        WHERE type='table' AND name=?
    """, (table,)).fetchone()

    if row is None or row[0] is None:
        raise RuntimeError(f"无法读取表结构: {table}")

    return row[0]


def get_index_sqls(conn: sqlite3.Connection, table: str):
    rows = conn.execute("""
        SELECT name, sql
        FROM sqlite_master
        WHERE type='index'
          AND tbl_name=?
          AND sql IS NOT NULL
        ORDER BY name
    """, (table,)).fetchall()

    return [r[1] for r in rows if r[1]]


def create_table_only(src_conn, dst_conn, table: str):
    create_sql = get_create_table_sql(src_conn, table)
    dst_conn.execute(create_sql)


def create_indexes(src_conn, dst_conn, table: str):
    index_sqls = get_index_sqls(src_conn, table)

    for sql in index_sqls:
        try:
            dst_conn.execute(sql)
        except sqlite3.OperationalError as e:
            log(f"  [警告] 创建索引失败，已跳过: {e}")
            log(f"         SQL: {sql}")


def vacuum_db(path: Path):
    conn = sqlite3.connect(path)
    try:
        conn.execute("VACUUM")
        conn.execute("ANALYZE")
    finally:
        conn.close()


def detect_gene_column(columns):
    """
    尽量自动识别哪个字段是基因 ID。
    兼容你之前数据库里出现过的字段名。
    """
    candidates = [
        "gene_id",
        "primary_gene_id",
        "gene",
        "gene.name",
        "cs_gene_id",
        "fielder_gene_id",
        "self_homolog_gene_id",
        "homolog_gene_id",
        "transcript_id",
    ]

    for c in candidates:
        if c in columns:
            return c

    # 兜底：找包含 gene 和 id 的列
    for c in columns:
        lc = c.lower()
        if "gene" in lc and "id" in lc:
            return c

    # 再兜底：找包含 gene 的列
    for c in columns:
        if "gene" in c.lower():
            return c

    return None


def parse_chr_from_gene_id(gene_id):
    """
    从小麦基因 ID 中解析染色体。

    支持：
    TraesCS1A02G000100
    TraesCS5B02G233300.1
    TraesFLD5B01G105200
    """
    if gene_id is None:
        return "unknown"

    s = str(gene_id)

    m = re.search(r"TraesCS([1-7][ABD])", s)
    if m:
        return m.group(1)

    m = re.search(r"TraesFLD([1-7][ABD])", s)
    if m:
        return m.group(1)

    return "unknown"


def copy_whole_table(src_conn, table: str, out_path: Path):
    """
    把一个表完整复制到一个新的 SQLite db。
    """
    if out_path.exists():
        out_path.unlink()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    log(f"  输出: {relpath(out_path)}")

    dst_conn = connect_sqlite(out_path)

    try:
        create_table_only(src_conn, dst_conn, table)

        columns = get_columns(src_conn, table)
        col_sql = ", ".join(quote_ident(c) for c in columns)
        placeholders = ", ".join(["?"] * len(columns))

        select_sql = f"SELECT {col_sql} FROM {quote_ident(table)}"
        insert_sql = f"INSERT INTO {quote_ident(table)} ({col_sql}) VALUES ({placeholders})"

        src_cur = src_conn.execute(select_sql)

        total = 0
        while True:
            rows = src_cur.fetchmany(BATCH_SIZE)
            if not rows:
                break

            dst_conn.executemany(insert_sql, rows)
            total += len(rows)

            if total % 50000 == 0:
                log(f"    已复制 {total} 行")

        dst_conn.commit()

        log("    创建索引...")
        create_indexes(src_conn, dst_conn, table)
        dst_conn.commit()

    finally:
        dst_conn.close()

    vacuum_db(out_path)

    size = file_size_mb(out_path)
    log(f"  完成: {table}, 行数 {total}, 大小 {size:.2f} MB")

    return {
        "strategy": "whole_table",
        "path": relpath(out_path),
        "rows": total,
        "size_mb": round(size, 2),
    }


def split_table_by_chr(src_conn, table: str, out_subdir: Path):
    """
    按染色体拆分表。
    每个染色体一个 SQLite db。
    """
    out_subdir.mkdir(parents=True, exist_ok=True)

    columns = get_columns(src_conn, table)
    gene_col = detect_gene_column(columns)

    if gene_col is None:
        raise RuntimeError(f"表 {table} 无法识别 gene_id 字段，不能按染色体拆分。字段为: {columns}")

    log(f"  识别到基因列: {gene_col}")

    col_sql = ", ".join(quote_ident(c) for c in columns)
    placeholders = ", ".join(["?"] * len(columns))

    select_sql = f"SELECT {col_sql} FROM {quote_ident(table)}"
    insert_sql = f"INSERT INTO {quote_ident(table)} ({col_sql}) VALUES ({placeholders})"

    gene_idx = columns.index(gene_col)

    shard_conns = {}
    shard_counts = {}
    shard_buffers = {}

    def ensure_shard(chr_name):
        if chr_name in shard_conns:
            return

        out_path = out_subdir / f"{table}_{chr_name}.db"
        if out_path.exists():
            out_path.unlink()

        log(f"  创建分库: {relpath(out_path)}")

        conn = connect_sqlite(out_path)
        create_table_only(src_conn, conn, table)

        shard_conns[chr_name] = {
            "conn": conn,
            "path": out_path,
        }
        shard_counts[chr_name] = 0
        shard_buffers[chr_name] = []

    def flush_shard(chr_name):
        buf = shard_buffers[chr_name]
        if not buf:
            return

        conn = shard_conns[chr_name]["conn"]
        conn.executemany(insert_sql, buf)
        shard_counts[chr_name] += len(buf)
        shard_buffers[chr_name] = []

    src_cur = src_conn.execute(select_sql)

    total = 0
    while True:
        rows = src_cur.fetchmany(BATCH_SIZE)
        if not rows:
            break

        for row in rows:
            gene_id = row[gene_idx]
            chr_name = parse_chr_from_gene_id(gene_id)

            ensure_shard(chr_name)
            shard_buffers[chr_name].append(row)

            if len(shard_buffers[chr_name]) >= BATCH_SIZE:
                flush_shard(chr_name)

        total += len(rows)

        if total % 50000 == 0:
            log(f"    已扫描 {total} 行")

    # flush all
    for chr_name in list(shard_buffers.keys()):
        flush_shard(chr_name)

    shards = {}

    # 创建索引、关闭、压缩
    for chr_name, info in shard_conns.items():
        conn = info["conn"]
        path = info["path"]

        log(f"  分库 {chr_name}: 创建索引...")
        create_indexes(src_conn, conn, table)
        conn.commit()
        conn.close()

        vacuum_db(path)

        size = file_size_mb(path)
        rows = shard_counts[chr_name]

        log(f"  完成分库: {path.name}, 行数 {rows}, 大小 {size:.2f} MB")

        shards[chr_name] = {
            "path": relpath(path),
            "rows": rows,
            "size_mb": round(size, 2),
        }

    return {
        "strategy": "by_chr",
        "gene_column": gene_col,
        "total_rows": total,
        "shards": shards,
    }


def check_large_files(out_dir: Path, max_mb: int = MAX_FILE_MB):
    large = []

    for p in out_dir.rglob("*.db"):
        size = file_size_mb(p)
        if size > max_mb:
            large.append((p, size))

    return large


def main():
    start = time.time()

    if not SOURCE_DB.exists():
        raise FileNotFoundError(f"找不到源数据库: {SOURCE_DB}")

    log("=" * 80)
    log("开始拆分 wheat_toolkit.db")
    log(f"项目目录: {PROJECT_ROOT}")
    log(f"源数据库: {SOURCE_DB}")
    log(f"源数据库大小: {file_size_mb(SOURCE_DB):.2f} MB")
    log(f"输出目录: {OUT_DIR}")
    log("=" * 80)

    if OUT_DIR.exists():
        log(f"检测到旧目录，先删除: {OUT_DIR}")
        shutil.rmtree(OUT_DIR)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    src_conn = sqlite3.connect(SOURCE_DB)

    try:
        tables = get_tables(src_conn)

        log("源数据库表：")
        for t in tables:
            log(f"  - {t}")

        manifest = {
            "source_db": "wheat_toolkit.db",
            "source_size_mb": round(file_size_mb(SOURCE_DB), 2),
            "max_file_mb": MAX_FILE_MB,
            "batch_size": BATCH_SIZE,
            "tables": {},
        }

        for table in tables:
            log("\n" + "=" * 80)
            log(f"处理表: {table}")
            log("=" * 80)

            if table in SPLIT_BY_CHR_TABLES:
                out_subdir = OUT_DIR / table
                info = split_table_by_chr(src_conn, table, out_subdir)

            elif table in CORE_TABLES:
                out_path = OUT_DIR / "core" / f"{table}.db"
                info = copy_whole_table(src_conn, table, out_path)

            elif table in HOMOLOG_TABLES:
                out_path = OUT_DIR / "homolog" / f"{table}.db"
                info = copy_whole_table(src_conn, table, out_path)

            else:
                out_path = OUT_DIR / "misc" / f"{table}.db"
                info = copy_whole_table(src_conn, table, out_path)

            manifest["tables"][table] = info

        manifest_path = OUT_DIR / "manifest.json"
        with manifest_path.open("w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)

    finally:
        src_conn.close()

    log("\n" + "=" * 80)
    log("拆库完成")
    log("=" * 80)

    log(f"manifest 文件: {relpath(OUT_DIR / 'manifest.json')}")

    log("\n生成的 db 文件：")
    for p in sorted(OUT_DIR.rglob("*.db")):
        log(f"  {relpath(p)}\t{file_size_mb(p):.2f} MB")

    large_files = check_large_files(OUT_DIR, MAX_FILE_MB)

    if large_files:
        log("\n[警告] 以下文件超过 95 MB，可能不适合上传 GitHub：")
        for p, size in large_files:
            log(f"  {relpath(p)}\t{size:.2f} MB")
        log("\n这些文件需要进一步二次拆分。把这部分输出发我，我继续帮你拆。")
    else:
        log("\n检查通过：所有 .db 文件均小于 95 MB。")

    elapsed = time.time() - start
    log(f"\n总耗时: {elapsed:.1f} 秒")


if __name__ == "__main__":
    main()
