import bisect
import json
import math
import random
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
    "p_value",
    "q_value",
    "significant",
    "p_level",
    "significance_rank",
    "confidence_level",
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


def load_precomputed_background(json_path):
    """
    读取离线预计算的 JASPAR background score distribution JSON。

    文件不存在或格式异常时返回 None，避免 Streamlit 页面因为缺少预计算
    文件直接崩溃。返回结构包含 background、n_samples 和 distributions。
    """
    path = Path(json_path)
    if not path.exists():
        return None

    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None

    distributions = payload.get("distributions")
    if not isinstance(distributions, dict):
        return None

    return {
        "name": payload.get("name", ""),
        "background_name": payload.get("background_name", ""),
        "background": payload.get("background") or {},
        "n_samples": payload.get("n_samples"),
        "motif_count": payload.get("motif_count"),
        "distributions": distributions,
    }


def load_precomputed_thresholds(json_path):
    """
    读取轻量级 JASPAR background threshold JSON。

    文件不存在或格式异常时返回 None，避免 Streamlit 因缺少阈值表而崩溃。
    """
    path = Path(json_path)
    if not path.exists():
        return None

    try:
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None

    thresholds = payload.get("thresholds")
    if not isinstance(thresholds, dict):
        return None

    return {
        "name": payload.get("name", ""),
        "background_name": payload.get("background_name", ""),
        "background": payload.get("background") or {},
        "n_samples": payload.get("n_samples"),
        "p_levels": payload.get("p_levels") or [],
        "motif_count": payload.get("motif_count"),
        "thresholds": thresholds,
    }


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


def estimate_background_from_records(records):
    """
    根据输入序列估计 A/C/G/T 背景频率。

    N 和其他非 A/C/G/T 字符会被忽略。如果有效碱基太少，则使用均匀背景，
    避免极端输入导致随机背景模型不稳定。
    """
    counts = {base: 0 for base in BASES}
    total = 0

    for sequence in (records or {}).values():
        for base in clean_sequence(sequence):
            if base in counts:
                counts[base] += 1
                total += 1

    if total < 20:
        return {base: 0.25 for base in BASES}

    return {base: counts[base] / total for base in BASES}


def generate_random_score_distribution(motif, background, n_samples=50000, seed=123):
    """
    使用给定背景模型随机生成 DNA window，并计算该 motif 的 PWM raw score 分布。

    返回升序排序后的分数列表。使用 random.Random(seed) 保证结果可复现。
    """
    try:
        motif_length = int(motif.get("length") or 0)
    except (TypeError, ValueError):
        return []

    pwm = motif.get("pwm") or {}
    if motif_length <= 0 or not pwm:
        return []

    bases = list(BASES)
    weights = []
    for base in bases:
        try:
            weights.append(float((background or {}).get(base, 0.25)))
        except (TypeError, ValueError):
            weights.append(0.25)

    weight_sum = sum(weights)
    if weight_sum <= 0:
        weights = [0.25, 0.25, 0.25, 0.25]
    else:
        weights = [weight / weight_sum for weight in weights]

    rng = random.Random(seed)
    scores = []
    for _ in range(int(n_samples)):
        # 按背景频率生成随机窗口，再复用现有 PWM 打分函数。
        window = "".join(rng.choices(bases, weights=weights, k=motif_length))
        score = score_window(window, pwm)
        if score is not None:
            scores.append(score)

    scores.sort()
    return scores


def pvalue_from_score_distribution(raw_score, sorted_scores):
    """
    根据随机背景分数分布计算右尾 p-value。

    p-value 表示背景模型中 score >= raw_score 的概率。加入 +1 平滑，
    避免模拟次数有限时出现 p-value 为 0。
    """
    if raw_score is None or not sorted_scores:
        return 1.0

    n = len(sorted_scores)
    index = bisect.bisect_left(sorted_scores, float(raw_score))
    n_ge = n - index
    return (n_ge + 1) / (n + 1)


def bh_fdr_correction(p_values):
    """
    Benjamini-Hochberg FDR 校正。

    输入 p-value 列表，返回同长度 q-value 列表。这里手动实现，
    不依赖 statsmodels，并通过从大到小回填保证 q-value 单调性。
    """
    n = len(p_values)
    if n == 0:
        return []

    cleaned = []
    for idx, p_value in enumerate(p_values):
        try:
            p = float(p_value)
            if math.isnan(p) or p < 0:
                p = 1.0
            p = min(p, 1.0)
        except (TypeError, ValueError):
            p = 1.0
        cleaned.append((idx, p))

    order = sorted(cleaned, key=lambda item: item[1])
    q_values = [1.0] * n
    running_min = 1.0

    for rank in range(n, 0, -1):
        original_idx, p = order[rank - 1]
        q = min(p * n / rank, 1.0)
        running_min = min(running_min, q)
        q_values[original_idx] = running_min

    return q_values


def assign_p_level(raw_score, motif_thresholds):
    """
    根据预计算 score cutoff 给 raw_score 分配离散 p-level。
    """
    if raw_score is None or not isinstance(motif_thresholds, dict):
        return ">0.05"

    ordered_levels = [0.0001, 0.001, 0.01, 0.05]
    for level in ordered_levels:
        cutoff = motif_thresholds.get(str(level))
        if cutoff is None:
            cutoff = motif_thresholds.get(f"{level:g}")
        try:
            if float(raw_score) >= float(cutoff):
                return f"<={level:g}"
        except (TypeError, ValueError):
            continue

    return ">0.05"


def p_level_to_rank(p_level: str) -> int:
    """
    将 p_level 字符串转换为排序等级，数字越小越显著。
    """
    order = {
        "<=0.0001": 1,
        "<=0.001": 2,
        "<=0.01": 3,
        "<=0.05": 4,
        ">0.05": 5,
    }
    return order.get(str(p_level), 5)


def p_level_passes_cutoff(p_level: str, cutoff: float = 0.01) -> bool:
    """
    判断离散 p_level 是否通过用户选择的 p-level cutoff。
    """
    cutoff_to_rank = {
        0.0001: 1,
        0.001: 2,
        0.01: 3,
        0.05: 4,
    }
    try:
        cutoff_value = float(cutoff)
    except (TypeError, ValueError):
        cutoff_value = 0.01
    target_rank = cutoff_to_rank.get(cutoff_value, 3)
    return p_level_to_rank(p_level) <= target_rank


def confidence_level_from_hit(p_level: str, rel_score) -> str:
    """
    根据 p-level 和 relative_score 给候选 TFBS 分级。
    """
    try:
        rel = float(rel_score)
    except (TypeError, ValueError):
        rel = 0.0

    rank = p_level_to_rank(p_level)
    if rank <= p_level_to_rank("<=0.001") and rel >= 0.95:
        return "high"
    if rank <= p_level_to_rank("<=0.01") and rel >= 0.90:
        return "medium"
    if rank <= p_level_to_rank("<=0.05"):
        return "low"
    return "not_significant"


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
        "p_value": None,
        "q_value": None,
        "significant": None,
        "p_level": ">0.05",
        "significance_rank": 5,
        "confidence_level": "not_significant",
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


def _stable_motif_seed(matrix_id: str, base_seed: int) -> int:
    """
    为每个 motif 生成稳定随机种子，避免不同 motif 共用完全相同的随机窗口序列。
    """
    offset = sum((idx + 1) * ord(char) for idx, char in enumerate(str(matrix_id)))
    return int(base_seed) + offset


def scan_sequences_with_jaspar_significance(
    records,
    motifs,
    relative_cutoff=0.85,
    qvalue_cutoff=0.05,
    scan_reverse=True,
    selected_matrix_ids=None,
    max_total_hits=20000,
    background_mode="input",
    n_background_samples=50000,
    random_seed=123,
    precomputed_distributions=None,
) -> pd.DataFrame:
    """
    在原有 PWM relative score 扫描基础上增加统计显著性评估。

    流程：
    1. 先用 relative_cutoff 做初筛，得到候选 hit。
    2. 根据输入序列或均匀模型建立背景碱基频率。
    3. 对每个候选 motif 模拟随机背景分数分布并计算 p-value。
    4. 对所有候选 hit 做 BH-FDR 校正，得到 q-value。
    """
    motifs = list(motifs)
    candidate_df = scan_sequences_with_jaspar(
        records=records,
        motifs=motifs,
        cutoff=relative_cutoff,
        scan_reverse=scan_reverse,
        selected_matrix_ids=selected_matrix_ids,
        max_total_hits=max_total_hits,
    )

    if candidate_df.empty:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    candidate_df = candidate_df.copy()

    precomputed_background = {}
    if isinstance(precomputed_distributions, dict) and "distributions" in precomputed_distributions:
        precomputed_background = precomputed_distributions.get("background") or {}
        precomputed_distributions = precomputed_distributions.get("distributions") or {}

    use_input_background = background_mode == "input"
    if use_input_background:
        background = estimate_background_from_records(records)
    elif precomputed_background:
        background = {base: float(precomputed_background.get(base, 0.25)) for base in BASES}
    elif background_mode == "uniform":
        background = {base: 0.25 for base in BASES}
    else:
        # 如果 CS/Fielder 预计算分布不存在且没有经验背景信息，则用均匀背景实时模拟作为兜底。
        background = {base: 0.25 for base in BASES}

    motif_lookup = {str(motif.get("matrix_id", "")): motif for motif in motifs}
    distribution_cache = {}
    p_values = []

    for matrix_id, group in candidate_df.groupby("matrix_id", sort=False):
        motif = motif_lookup.get(str(matrix_id))
        if motif is None:
            p_values.extend([1.0] * len(group))
            continue

        if matrix_id not in distribution_cache:
            # 优先使用离线预计算分布；只有缺失或 input background 模式才实时模拟。
            precomputed_scores = None
            if not use_input_background and isinstance(precomputed_distributions, dict):
                precomputed_scores = precomputed_distributions.get(str(matrix_id))

            if precomputed_scores:
                distribution_cache[matrix_id] = precomputed_scores
            else:
                distribution_cache[matrix_id] = generate_random_score_distribution(
                    motif=motif,
                    background=background,
                    n_samples=n_background_samples,
                    seed=_stable_motif_seed(matrix_id, random_seed),
                )

        sorted_scores = distribution_cache[matrix_id]
        for raw_score_value in group["raw_score"].tolist():
            p_values.append(pvalue_from_score_distribution(raw_score_value, sorted_scores))

    # groupby 遍历会按分组顺序收集 p-value，需要按同样顺序写回原始索引。
    ordered_indices = []
    for _, group in candidate_df.groupby("matrix_id", sort=False):
        ordered_indices.extend(group.index.tolist())

    p_value_by_index = dict(zip(ordered_indices, p_values))
    candidate_df["p_value"] = [p_value_by_index.get(idx, 1.0) for idx in candidate_df.index]
    candidate_df["q_value"] = bh_fdr_correction(candidate_df["p_value"].tolist())
    candidate_df["significant"] = candidate_df["q_value"] <= float(qvalue_cutoff)

    return candidate_df[RESULT_COLUMNS]


def scan_sequences_with_jaspar_thresholds(
    records,
    motifs,
    relative_cutoff=0.90,
    p_level_cutoff=0.01,
    scan_reverse=True,
    selected_matrix_ids=None,
    max_total_hits=20000,
    precomputed_thresholds=None,
    top_n_per_motif_sequence=3,
) -> pd.DataFrame:
    """
    使用轻量级预计算 score cutoff 表对 PWM hits 做显著性分级。

    该函数不加载完整背景分布，也不计算精确 p-value/q-value，适合作为
    Streamlit 默认快速模式。
    """
    motifs = list(motifs)
    candidate_df = scan_sequences_with_jaspar(
        records=records,
        motifs=motifs,
        cutoff=relative_cutoff,
        scan_reverse=scan_reverse,
        selected_matrix_ids=selected_matrix_ids,
        max_total_hits=max_total_hits,
    )

    if candidate_df.empty:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    candidate_df = candidate_df.copy()
    if isinstance(precomputed_thresholds, dict) and "thresholds" in precomputed_thresholds:
        threshold_map = precomputed_thresholds.get("thresholds") or {}
    elif isinstance(precomputed_thresholds, dict):
        threshold_map = precomputed_thresholds
    else:
        threshold_map = {}

    p_levels = []
    ranks = []
    confidence_levels = []
    significant_flags = []

    for _, row in candidate_df.iterrows():
        matrix_id = str(row.get("matrix_id", ""))
        motif_thresholds = threshold_map.get(matrix_id, {})
        p_level = assign_p_level(row.get("raw_score"), motif_thresholds)
        rank = p_level_to_rank(p_level)
        confidence = confidence_level_from_hit(p_level, row.get("relative_score"))

        p_levels.append(p_level)
        ranks.append(rank)
        confidence_levels.append(confidence)
        significant_flags.append(p_level_passes_cutoff(p_level, p_level_cutoff))

    candidate_df["p_level"] = p_levels
    candidate_df["significance_rank"] = ranks
    candidate_df["confidence_level"] = confidence_levels
    candidate_df["significant"] = significant_flags
    candidate_df["p_value"] = None
    candidate_df["q_value"] = None

    candidate_df = candidate_df.sort_values(
        ["sequence_id", "matrix_id", "significance_rank", "relative_score", "raw_score"],
        ascending=[True, True, True, False, False],
    )

    if top_n_per_motif_sequence is not None and int(top_n_per_motif_sequence) > 0:
        candidate_df = (
            candidate_df.groupby(["sequence_id", "matrix_id"], dropna=False)
            .head(int(top_n_per_motif_sequence))
            .reset_index(drop=True)
        )

    candidate_df = candidate_df.sort_values(
        ["significance_rank", "relative_score", "raw_score"],
        ascending=[True, False, False],
    ).reset_index(drop=True)

    return candidate_df[RESULT_COLUMNS]
