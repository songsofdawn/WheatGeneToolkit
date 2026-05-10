import json
import sqlite3
from pathlib import Path


BASES = ("A", "C", "G", "T")


def count_bases(sequence: str, counts: dict):
    """
    统计一条 promoter 序列中的 A/C/G/T，忽略 N 和其他字符。
    """
    total = 0
    for base in (sequence or "").upper():
        if base in counts:
            counts[base] += 1
            total += 1
    return total


def find_sqlite_files(project_root: Path, directory_name: str, single_db_name: str):
    """
    自动识别拆分目录或单个 sqlite 数据库。
    """
    data_dir = project_root / "data"
    candidates = []

    split_dir = data_dir / "db" / directory_name
    if split_dir.exists():
        candidates.extend(sorted(split_dir.glob("*.db")))

    for path in [
        data_dir / "db" / single_db_name,
        data_dir / single_db_name,
        data_dir / directory_name / single_db_name,
    ]:
        if path.exists():
            candidates.append(path)

    return sorted(set(candidates))


def find_sequence_table_and_column(db_path: Path):
    """
    在 sqlite 数据库中查找包含 promoter_sequence 列的表。
    """
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [row[0] for row in cursor.fetchall() if not row[0].startswith("sqlite_")]

        for table in tables:
            cursor.execute(f"PRAGMA table_info({table})")
            columns = [row[1] for row in cursor.fetchall()]
            lower_columns = {col.lower(): col for col in columns}
            for candidate in ["promoter_sequence", "sequence", "seq"]:
                if candidate in lower_columns:
                    return table, lower_columns[candidate]

    return None, None


def compute_background_from_sqlite_files(db_files):
    """
    从一组 promoter sqlite 数据库中统计经验背景频率。
    """
    counts = {base: 0 for base in BASES}
    total_bases = 0
    sequence_count = 0
    used_files = []

    for db_path in db_files:
        table, sequence_col = find_sequence_table_and_column(db_path)
        if table is None:
            print(f"[WARN] 未在数据库中找到 promoter 序列表，跳过: {db_path}")
            continue

        print(f"[INFO] 读取 {db_path} :: {table}.{sequence_col}")
        used_files.append(str(db_path))
        with sqlite3.connect(db_path) as conn:
            cursor = conn.cursor()
            cursor.execute(f"SELECT {sequence_col} FROM {table}")
            for (sequence,) in cursor:
                total_bases += count_bases(sequence, counts)
                sequence_count += 1

    if total_bases == 0:
        return None

    background = {base: counts[base] / total_bases for base in BASES}
    return {
        "counts": counts,
        "background": background,
        "sequence_count": sequence_count,
        "total_bases": total_bases,
        "used_files": used_files,
    }


def write_background(output_path: Path, name: str, genome: str, result: dict):
    """
    写出 promoter empirical background JSON。
    """
    payload = {
        "name": name,
        "source": "computed from local promoter sequence database",
        "genome": genome,
        "promoter_length": 2000,
        "sequence_count": result["sequence_count"],
        "total_bases": result["total_bases"],
        "counts": result["counts"],
        "background": result["background"],
        "used_files": result["used_files"],
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    print(f"[OK] 已输出: {output_path}")
    print(f"[OK] background: {payload['background']}")


def build_one(project_root: Path, label: str, directory_name: str, single_db_name: str, output_name: str, genome: str):
    """
    构建一个基因组的 promoter empirical background。
    """
    print("=" * 80)
    print(f"开始识别 {label} promoter 数据库")
    db_files = find_sqlite_files(project_root, directory_name, single_db_name)
    if not db_files:
        print(f"[INFO] 未找到 {label} promoter 数据库，已跳过。候选路径包括:")
        print(f"  data/db/{directory_name}/")
        print(f"  data/db/{single_db_name}")
        print(f"  data/{single_db_name}")
        return

    print(f"[INFO] 找到 {len(db_files)} 个数据库文件")
    result = compute_background_from_sqlite_files(db_files)
    if result is None:
        print(f"[WARN] {label} promoter 数据库没有可统计的 A/C/G/T 碱基，已跳过。")
        return

    result["used_files"] = [
        str(Path(path).relative_to(project_root)).replace("\\", "/")
        for path in result["used_files"]
    ]

    jaspar_dir = project_root / "data" / "motif_db" / "jaspar_plants"
    write_background(
        output_path=jaspar_dir / output_name,
        name=f"{label} promoter 2000bp empirical background",
        genome=genome,
        result=result,
    )


def main():
    project_root = Path(__file__).resolve().parents[1]
    build_one(
        project_root=project_root,
        label="Chinese Spring",
        directory_name="gene_promoter_sequence",
        single_db_name="gene_promoter_sequence.db",
        output_name="background_cs_promoter_2000.json",
        genome="Chinese Spring",
    )
    build_one(
        project_root=project_root,
        label="Fielder",
        directory_name="fielder_promoter_sequence",
        single_db_name="fielder_promoter_sequence.db",
        output_name="background_fielder_promoter_2000.json",
        genome="Fielder",
    )


if __name__ == "__main__":
    main()
