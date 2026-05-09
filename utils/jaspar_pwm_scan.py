import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd


BASES = ("A", "C", "G", "T")
RESULT_COLUMNS = [
    "sequence_id",
    "matrix_id",
    "tf_name",
    "consensus",
    "motif_length",
    "start",
    "end",
    "strand",
    "matched_seq",
    "raw_score",
    "relative_score",
    "distance_to_sequence_end",
    "species",
    "tax_group",
    "family",
    "class",
    "collection",
]


def clean_sequence(seq: str) -> str:
    """
    清理 DNA 序列，只保留 A/C/G/T/N，并统一转成大写。
    """
    if not seq:
        return ""
    allowed = set("ACGTN")
    return "".join(base for base in seq.upper() if base in allowed)


def reverse_complement(seq: str) -> str:
    """
    返回 DNA 序列的反向互补序列。
    """
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1].upper()


def parse_fasta_or_plain(text: str) -> Dict[str, str]:
    """
    解析 FASTA 或纯 DNA 输入。

    - FASTA: 返回 {序列名: 清理后的序列}
    - 纯序列: 返回 {"input_sequence_1": 清理后的序列}
    """
    text = (text or "").strip()
    if not text:
        return {}

    lines = [line.strip() for line in text.splitlines() if line.strip()]
    is_fasta = any(line.startswith(">") for line in lines)

    if not is_fasta:
        return {"input_sequence_1": clean_sequence(text)}

    records: Dict[str, str] = {}
    current_name: Optional[str] = None
    current_chunks: List[str] = []

    def flush_record():
        nonlocal current_name, current_chunks
        if current_name is None:
            return
        seq = clean_sequence("".join(current_chunks))
        if seq:
            name = current_name
            # 避免重复 FASTA header 覆盖前面的序列。
            if name in records:
                suffix = 2
                while f"{name}_{suffix}" in records:
                    suffix += 1
                name = f"{name}_{suffix}"
            records[name] = seq
        current_name = None
        current_chunks = []

    for line in lines:
        if line.startswith(">"):
            flush_record()
            header = line[1:].strip()
            current_name = header or f"input_sequence_{len(records) + 1}"
            current_chunks = []
        elif current_name is not None:
            current_chunks.append(line)

    flush_record()
    return records


def load_jaspar_pwm(json_path) -> List[dict]:
    """
    读取 build_jaspar_pwm.py 生成的 jaspar_plants_pwm.json。
    """
    path = Path(json_path)
    with path.open("r", encoding="utf-8") as handle:
        motifs = json.load(handle)
    if not isinstance(motifs, list):
        raise ValueError(f"JASPAR PWM JSON 格式异常，预期为 motif list: {path}")
    return motifs


def score_window(window: str, pwm: dict):
    """
    计算一个窗口的 PWM raw score。

    如果窗口包含 N 或非 A/C/G/T，则返回 None，避免给不确定碱基打分。
    """
    window = (window or "").upper()
    if any(base not in BASES for base in window):
        return None

    score = 0.0
    for idx, base in enumerate(window):
        try:
            score += float(pwm[base][idx])
        except (KeyError, IndexError, TypeError, ValueError):
            return None
    return score


def relative_score(raw_score, min_score, max_score):
    """
    计算 JASPAR 常用的 relative score:
    (raw_score - min_score) / (max_score - min_score)
    """
    try:
        denominator = float(max_score) - float(min_score)
        if denominator == 0:
            return None
        return (float(raw_score) - float(min_score)) / denominator
    except (TypeError, ValueError):
        return None


def _metadata_value(motif: dict, key: str) -> str:
    metadata = motif.get("metadata") or {}
    return str(metadata.get(key, "") or "")


def _motif_name(motif: dict) -> str:
    return str(motif.get("name") or _metadata_value(motif, "name") or "")


def _make_hit_row(
    sequence_id: str,
    motif: dict,
    seq_len: int,
    start: int,
    end: int,
    strand: str,
    matched_seq: str,
    raw_score_value: float,
    relative_score_value: float,
) -> dict:
    """
    组装输出行。坐标使用 1-based，并保留常用 metadata 字段。
    """
    return {
        "sequence_id": sequence_id,
        "matrix_id": motif.get("matrix_id", ""),
        "tf_name": _motif_name(motif),
        "consensus": motif.get("consensus", ""),
        "motif_length": int(motif.get("length") or len(matched_seq)),
        "start": start,
        "end": end,
        "strand": strand,
        "matched_seq": matched_seq,
        "raw_score": raw_score_value,
        "relative_score": relative_score_value,
        # 假设输入为 ATG 上游 promoter，则 start - seq_len - 1 近似为相对 ATG 距离。
        "distance_to_sequence_end": start - seq_len - 1,
        "species": _metadata_value(motif, "species"),
        "tax_group": _metadata_value(motif, "tax_group"),
        "family": _metadata_value(motif, "family"),
        "class": _metadata_value(motif, "class"),
        "collection": _metadata_value(motif, "collection"),
    }


def scan_sequence_with_motif(
    sequence_id: str,
    sequence: str,
    motif: dict,
    cutoff: float = 0.85,
    scan_reverse: bool = True,
    max_hits: Optional[int] = None,
) -> List[dict]:
    """
    用单个 motif 扫描一条序列，支持正链和反向互补链。
    """
    sequence = clean_sequence(sequence)
    seq_len = len(sequence)
    motif_length = int(motif.get("length") or 0)
    pwm = motif.get("pwm") or {}
    min_score = motif.get("min_score")
    max_score = motif.get("max_score")
    hits: List[dict] = []

    if not sequence or motif_length <= 0 or seq_len < motif_length:
        return hits

    def add_hit(start: int, end: int, strand: str, matched_seq: str, raw_score_value: float, rel_value: float):
        hits.append(
            _make_hit_row(
                sequence_id=sequence_id,
                motif=motif,
                seq_len=seq_len,
                start=start,
                end=end,
                strand=strand,
                matched_seq=matched_seq,
                raw_score_value=raw_score_value,
                relative_score_value=rel_value,
            )
        )

    for idx in range(seq_len - motif_length + 1):
        window = sequence[idx : idx + motif_length]
        raw = score_window(window, pwm)
        if raw is None:
            continue
        rel = relative_score(raw, min_score, max_score)
        if rel is not None and rel >= cutoff:
            add_hit(
                start=idx + 1,
                end=idx + motif_length,
                strand="+",
                matched_seq=window,
                raw_score_value=raw,
                rel_value=rel,
            )
            if max_hits is not None and len(hits) >= max_hits:
                return hits

    if not scan_reverse:
        return hits

    rc_sequence = reverse_complement(sequence)
    for rc_idx in range(seq_len - motif_length + 1):
        window = rc_sequence[rc_idx : rc_idx + motif_length]
        raw = score_window(window, pwm)
        if raw is None:
            continue
        rel = relative_score(raw, min_score, max_score)
        if rel is not None and rel >= cutoff:
            # 反向互补链坐标换算回原始输入序列坐标，仍使用 1-based。
            start = seq_len - (rc_idx + motif_length) + 1
            end = seq_len - rc_idx
            add_hit(
                start=start,
                end=end,
                strand="-",
                matched_seq=reverse_complement(window),
                raw_score_value=raw,
                rel_value=rel,
            )
            if max_hits is not None and len(hits) >= max_hits:
                return hits

    return hits


def scan_sequences_with_jaspar(
    records: Dict[str, str],
    motifs: Iterable[dict],
    cutoff: float = 0.85,
    scan_reverse: bool = True,
    selected_matrix_ids=None,
    max_total_hits: int = 20000,
) -> pd.DataFrame:
    """
    用所有 JASPAR motifs 扫描多条序列，并返回 pandas.DataFrame。
    """
    if not records:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    selected_ids = set(selected_matrix_ids or [])
    use_filter = bool(selected_ids)
    all_hits: List[dict] = []

    for motif in motifs:
        matrix_id = str(motif.get("matrix_id", ""))
        if use_filter and matrix_id not in selected_ids:
            continue

        for sequence_id, sequence in records.items():
            remaining = max_total_hits - len(all_hits)
            if remaining <= 0:
                return pd.DataFrame(all_hits, columns=RESULT_COLUMNS)

            hits = scan_sequence_with_motif(
                sequence_id=sequence_id,
                sequence=sequence,
                motif=motif,
                cutoff=cutoff,
                scan_reverse=scan_reverse,
                max_hits=remaining,
            )
            all_hits.extend(hits)

            if len(all_hits) >= max_total_hits:
                return pd.DataFrame(all_hits, columns=RESULT_COLUMNS)

    return pd.DataFrame(all_hits, columns=RESULT_COLUMNS)
