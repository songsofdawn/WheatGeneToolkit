# -*- coding: utf-8 -*-

"""
二次拆分 gene_structure_feature 表。

原因：
第一次拆库后 data/db/core/gene_structure_feature.db 仍有 372 MB，
超过 GitHub 单文件限制，因此需要按染色体进一步拆分。

运行方式：
cd <project-root>
python scripts\\split_gene_structure_feature.py
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
OUT_SUBDIR = OUT_DIR / "gene_structure_feature"
OLD_BIG_DB = OUT_DIR / "core" / "gene_structure_feature.db"
MANIFEST_PATH = OUT_DIR / "manifest.json"

TABLE = "gene_structure_feature"
BATCH_SIZE = 10000
MAX_FILE_MB = 95


def log(msg):
    print(msg, flush=True)


def quote_ident(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def relpath(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT)).replace("\\", "/")


def file_size_mb(path: Path) -> float:
    if not path.exists():
        return 0.0
    return path.stat().st_size / 1024 / 1024


def connect_sqlite(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=MEMORY")
    conn.execute("PRAGMA locking_mode=EXCLUSIVE")
    return conn


def get_columns(conn, table):
    rows = conn.execute(f"PRAGMA table_info({quote_ident(table)})").fetchall()
    return [r[1] for r in rows]


def get_create_table_sql(conn, table):
    row = conn.execute("""
        SELECT sql
        FROM sqlite_master
        WHERE type='table' AND name=?
    """, (table,)).fetchone()

    if row is None or row[0] is None:
        raise RuntimeError(f"无法读取表结构: {table}")

    return row[0]


def get_index_sqls(conn, table):
    rows = conn.execute("""
        SELECT sql
        FROM sqlite_master
        WHERE type='index'
          AND tbl_name=?
          AND sql IS NOT NULL
        ORDER BY name
    """, (table,)).fetchall()

    return [r[0] for r in rows if r[0]]


def create_table_only(src_conn, dst_conn, table):
    dst_conn.execute(get_create_table_sql(src_conn, table))


def create_indexes(src_conn, dst_conn, table):
    for sql in get_index_sqls(src_conn, table):
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
    自动识别 gene_structure_feature 中的基因列。
    """
    candidates = [
        "primary_gene_id",
        "gene_id",
        "gene",
        "gene.name",
        "cs_gene_id",
        "fielder_gene_id",
        "transcript_id",
    ]

    for c in candidates:
        if c in columns:
            return c

    for c in columns:
        lc = c.lower()
        if "gene" in lc and "id" in lc:
            return c

    for c in columns:
        if "gene" in c.lower():
            return c

    raise RuntimeError(f"无法识别基因列，当前字段为: {columns}")


def parse_chr_from_gene_id(gene_id):
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


def update_manifest(gene_column, total_rows, shards):
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"找不到 manifest.json: {MANIFEST_PATH}")

    with MANIFEST_PATH.open("r", encoding="utf-8") as f:
        manifest = json.load(f)

    manifest["tables"][TABLE] = {
        "strategy": "by_chr",
        "gene_column": gene_column,
        "total_rows": total_rows,
        "shards": shards,
    }

    with MANIFEST_PATH.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)

    log(f"已更新 manifest: {relpath(MANIFEST_PATH)}")


def main():
    start = time.time()

    if not SOURCE_DB.exists():
        raise FileNotFoundError(f"找不到源数据库: {SOURCE_DB}")

    log("=" * 80)
    log("开始二次拆分 gene_structure_feature")
    log(f"源数据库: {SOURCE_DB}")
    log(f"输出目录: {OUT_SUBDIR}")
    log("=" * 80)

    if OUT_SUBDIR.exists():
        log(f"检测到旧目录，先删除: {relpath(OUT_SUBDIR)}")
        shutil.rmtree(OUT_SUBDIR)

    OUT_SUBDIR.mkdir(parents=True, exist_ok=True)

    src_conn = sqlite3.connect(SOURCE_DB)

    try:
        columns = get_columns(src_conn, TABLE)
        gene_column = detect_gene_column(columns)

        log(f"识别到基因列: {gene_column}")
        log(f"字段列表: {columns}")

        col_sql = ", ".join(quote_ident(c) for c in columns)
        placeholders = ", ".join(["?"] * len(columns))

        select_sql = f"SELECT {col_sql} FROM {quote_ident(TABLE)}"
        insert_sql = f"INSERT INTO {quote_ident(TABLE)} ({col_sql}) VALUES ({placeholders})"

        gene_idx = columns.index(gene_column)

        shard_conns = {}
        shard_buffers = {}
        shard_counts = {}

        def ensure_shard(chr_name):
            if chr_name in shard_conns:
                return

            out_path = OUT_SUBDIR / f"{TABLE}_{chr_name}.db"

            log(f"创建分库: {relpath(out_path)}")

            dst_conn = connect_sqlite(out_path)
            create_table_only(src_conn, dst_conn, TABLE)

            shard_conns[chr_name] = {
                "conn": dst_conn,
                "path": out_path,
            }
            shard_buffers[chr_name] = []
            shard_counts[chr_name] = 0

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

            if total % 100000 == 0:
                log(f"已扫描 {total} 行")

        for chr_name in list(shard_buffers.keys()):
            flush_shard(chr_name)

        shards = {}

        for chr_name, info in shard_conns.items():
            conn = info["conn"]
            path = info["path"]

            log(f"分库 {chr_name}: 创建索引...")
            create_indexes(src_conn, conn, TABLE)
            conn.commit()
            conn.close()

            vacuum_db(path)

            size = file_size_mb(path)
            rows = shard_counts[chr_name]

            log(f"完成分库: {path.name}, 行数 {rows}, 大小 {size:.2f} MB")

            shards[chr_name] = {
                "path": relpath(path),
                "rows": rows,
                "size_mb": round(size, 2),
            }

        update_manifest(gene_column, total, shards)

    finally:
        src_conn.close()

    log("\n检查拆分后文件大小：")

    too_large = []

    for p in sorted(OUT_SUBDIR.glob("*.db")):
        size = file_size_mb(p)
        log(f"  {relpath(p)}\t{size:.2f} MB")

        if size > MAX_FILE_MB:
            too_large.append((p, size))

    if too_large:
        log("\n[警告] 以下文件仍然超过 95 MB：")
        for p, size in too_large:
            log(f"  {relpath(p)}\t{size:.2f} MB")
        log("需要继续按更细粒度拆分。")
    else:
        log("\n检查通过：gene_structure_feature 所有分库均小于 95 MB。")

        if OLD_BIG_DB.exists():
            log(f"\n删除旧的大文件: {relpath(OLD_BIG_DB)}")
            OLD_BIG_DB.unlink()

    elapsed = time.time() - start

    log("=" * 80)
    log("gene_structure_feature 二次拆分完成")
    log(f"总耗时: {elapsed:.1f} 秒")
    log("=" * 80)


if __name__ == "__main__":
    main()
