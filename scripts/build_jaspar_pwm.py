# scripts/build_jaspar_pwm.py

import json
import math
import re
from pathlib import Path

import pandas as pd


BASES = ["A", "C", "G", "T"]


def parse_jaspar_pfm(pfm_path: Path):
    """
    解析 JASPAR 格式 PFM 文件。

    格式一般类似：
    >MA0001.1 AGL3
    A [ 1 2 3 4 ]
    C [ 0 1 0 2 ]
    G [ 3 0 1 0 ]
    T [ 0 1 2 1 ]
    """
    motifs = []

    current_id = None
    current_name = None
    current_counts = {}

    def flush_current():
        nonlocal current_id, current_name, current_counts

        if current_id is None:
            return

        if all(base in current_counts for base in BASES):
            length_set = {len(current_counts[base]) for base in BASES}
            if len(length_set) != 1:
                raise ValueError(f"PFM length mismatch for {current_id}")

            motifs.append({
                "matrix_id": current_id,
                "name": current_name,
                "counts": current_counts,
                "length": list(length_set)[0],
            })

        current_id = None
        current_name = None
        current_counts = {}

    with pfm_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                flush_current()

                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                current_id = parts[0]
                current_name = parts[1] if len(parts) > 1 else parts[0]
                current_counts = {}
                continue

            m = re.match(r"^([ACGT])\s*\[\s*(.*?)\s*\]\s*$", line)
            if m:
                base = m.group(1)
                values = [float(x) for x in m.group(2).split()]
                current_counts[base] = values

    flush_current()

    return motifs


def find_column(df: pd.DataFrame, candidates):
    """
    在 metadata 表里尽量鲁棒地寻找列名。
    """
    lower_map = {col.lower(): col for col in df.columns}

    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]

    return None


def load_metadata(metadata_path: Path):
    """
    读取 JASPAR metadata。
    不强依赖具体列名，尽量兼容不同版本。
    """
    if not metadata_path.exists():
        return {}

    df = pd.read_csv(metadata_path, sep="\t", dtype=str)
    df = df.fillna("")

    matrix_col = find_column(df, [
        "matrix_id",
        "matrix id",
        "matrix",
        "jaspar_id",
        "id",
    ])

    if matrix_col is None:
        print("[WARN] 未找到 matrix_id 列，metadata 将不会合并。")
        print("[WARN] metadata columns:", list(df.columns))
        return {}

    metadata = {}

    for _, row in df.iterrows():
        matrix_id = str(row[matrix_col]).strip()
        if not matrix_id:
            continue

        metadata[matrix_id] = {col: str(row[col]) for col in df.columns}

    return metadata


def pfm_to_pwm(counts, pseudocount=0.8, background=None):
    """
    PFM counts 转 PWM。

    counts:
        {
            "A": [...],
            "C": [...],
            "G": [...],
            "T": [...]
        }

    PWM 公式：
        p_base = (count_base + pseudocount * bg_base) / (total + pseudocount)
        pwm = log2(p_base / bg_base)
    """
    if background is None:
        background = {
            "A": 0.25,
            "C": 0.25,
            "G": 0.25,
            "T": 0.25,
        }

    length = len(counts["A"])

    pwm = {base: [] for base in BASES}
    prob_matrix = {base: [] for base in BASES}

    for i in range(length):
        total = sum(counts[base][i] for base in BASES)

        for base in BASES:
            bg = background[base]
            prob = (counts[base][i] + pseudocount * bg) / (total + pseudocount)
            score = math.log2(prob / bg)

            prob_matrix[base].append(prob)
            pwm[base].append(score)

    max_score = 0.0
    min_score = 0.0

    for i in range(length):
        col_scores = [pwm[base][i] for base in BASES]
        max_score += max(col_scores)
        min_score += min(col_scores)

    consensus = ""

    for i in range(length):
        best_base = max(BASES, key=lambda b: prob_matrix[b][i])
        consensus += best_base

    return {
        "pwm": pwm,
        "probability_matrix": prob_matrix,
        "min_score": min_score,
        "max_score": max_score,
        "consensus": consensus,
    }


def main():
    project_root = Path(__file__).resolve().parents[1]

    jaspar_dir = project_root / "data" / "motif_db" / "jaspar_plants"

    pfm_path = jaspar_dir / "JASPAR2026_CORE_plants_non-redundant_pfms_jaspar.txt"
    metadata_path = jaspar_dir / "ultimate_metadata_table_CORE.tsv"
    output_path = jaspar_dir / "jaspar_plants_pwm.json"

    print("=" * 80)
    print("JASPAR Plants PFM to PWM")
    print("=" * 80)
    print(f"PFM file:      {pfm_path}")
    print(f"Metadata file: {metadata_path}")
    print(f"Output file:   {output_path}")
    print("=" * 80)

    if not pfm_path.exists():
        raise FileNotFoundError(f"找不到 PFM 文件: {pfm_path}")

    motifs = parse_jaspar_pfm(pfm_path)
    print(f"读取到 PFM motif 数量: {len(motifs)}")

    metadata = load_metadata(metadata_path)
    print(f"读取到 metadata 条目数: {len(metadata)}")

    output = []

    for motif in motifs:
        matrix_id = motif["matrix_id"]

        pwm_info = pfm_to_pwm(
            motif["counts"],
            pseudocount=0.8,
            background={
                "A": 0.25,
                "C": 0.25,
                "G": 0.25,
                "T": 0.25,
            },
        )

        item = {
            "matrix_id": matrix_id,
            "name": motif["name"],
            "length": motif["length"],
            "counts": motif["counts"],
            "pwm": pwm_info["pwm"],
            "probability_matrix": pwm_info["probability_matrix"],
            "min_score": pwm_info["min_score"],
            "max_score": pwm_info["max_score"],
            "consensus": pwm_info["consensus"],
            "metadata": metadata.get(matrix_id, {}),
        }

        output.append(item)

    with output_path.open("w", encoding="utf-8") as f:
        json.dump(output, f, ensure_ascii=False)

    print("=" * 80)
    print("转换完成")
    print(f"输出文件: {output_path}")
    print(f"输出 motif 数量: {len(output)}")
    print("=" * 80)

    if output:
        print("示例 motif:")
        example = output[0]
        print("matrix_id:", example["matrix_id"])
        print("name:", example["name"])
        print("length:", example["length"])
        print("consensus:", example["consensus"])
        print("min_score:", example["min_score"])
        print("max_score:", example["max_score"])


if __name__ == "__main__":
    main()