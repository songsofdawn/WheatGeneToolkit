import argparse
import json
import math
import random
import time
from pathlib import Path


BASES = ("A", "C", "G", "T")
DEFAULT_P_LEVELS = [0.05, 0.01, 0.001, 0.0001]


def score_window(window: str, pwm: dict):
    """
    计算随机窗口的 PWM raw score；窗口异常时返回 None。
    """
    score = 0.0
    for idx, base in enumerate(window):
        try:
            score += float(pwm[base][idx])
        except (KeyError, IndexError, TypeError, ValueError):
            return None
    return score


def stable_motif_seed(matrix_id: str, base_seed: int) -> int:
    """
    为每个 motif 生成稳定随机种子，保证预计算结果可复现。
    """
    return int(base_seed) + sum((idx + 1) * ord(char) for idx, char in enumerate(str(matrix_id)))


def simulate_scores(motif: dict, background: dict, n_samples: int, seed: int):
    """
    按背景频率随机生成窗口，并返回升序 PWM raw score 列表。

    注意：这个完整分布只在离线阶段短暂存在于内存中，不写入 JSON。
    """
    try:
        motif_length = int(motif.get("length") or 0)
    except (TypeError, ValueError):
        return []

    pwm = motif.get("pwm") or {}
    if motif_length <= 0 or not pwm:
        return []

    bases = list(BASES)
    weights = [float(background.get(base, 0.25)) for base in bases]
    weight_sum = sum(weights)
    if weight_sum <= 0:
        weights = [0.25, 0.25, 0.25, 0.25]
    else:
        weights = [value / weight_sum for value in weights]

    rng = random.Random(seed)
    scores = []
    for _ in range(n_samples):
        window = "".join(rng.choices(bases, weights=weights, k=motif_length))
        score = score_window(window, pwm)
        if score is not None:
            scores.append(score)

    scores.sort()
    return scores


def score_thresholds_from_distribution(sorted_scores, p_levels):
    """
    从升序背景分数中提取少量右尾 p-level 对应的 score cutoff。

    对于 p_level，取最高 tail_count 个背景窗口的最低分数作为 cutoff，
    使 P(background_score >= cutoff) 近似为 p_level。
    """
    if not sorted_scores:
        return {}

    n = len(sorted_scores)
    thresholds = {}
    for p_level in p_levels:
        tail_count = max(1, int(math.ceil(float(p_level) * n)))
        index = max(0, n - tail_count)
        thresholds[str(p_level)] = sorted_scores[index]
    return thresholds


def load_background_json(path: Path):
    """
    读取 promoter empirical background JSON 中的 background 字段。
    """
    if not path.exists():
        return None

    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    background = payload.get("background")
    if not isinstance(background, dict):
        print(f"[WARN] background 字段不存在或格式异常，跳过: {path}")
        return None

    return {base: float(background.get(base, 0.25)) for base in BASES}


def build_one_background(
    motifs,
    background_name: str,
    background: dict,
    output_path: Path,
    n_samples: int,
    base_seed: int,
    p_levels,
):
    """
    为一个背景模型预计算所有 motif 的 p-level score cutoff。
    """
    start_time = time.time()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print(f"Background: {background_name}")
    print(f"Output:     {output_path}")
    print(f"Samples:    {n_samples}")
    print(f"P-levels:   {p_levels}")
    print(f"Background: {background}")
    print("=" * 80)

    thresholds = {}
    motif_count = len(motifs)
    for idx, motif in enumerate(motifs, 1):
        matrix_id = str(motif.get("matrix_id", ""))
        print(f"[{idx}/{motif_count}] {matrix_id} {motif.get('name', '')}")
        scores = simulate_scores(
            motif=motif,
            background=background,
            n_samples=n_samples,
            seed=stable_motif_seed(matrix_id, base_seed),
        )
        thresholds[matrix_id] = score_thresholds_from_distribution(scores, p_levels)

    payload = {
        "name": f"JASPAR Plants {background_name} background score thresholds",
        "background_name": background_name,
        "background": background,
        "n_samples": n_samples,
        "p_levels": p_levels,
        "motif_count": motif_count,
        "thresholds": thresholds,
    }

    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    elapsed = time.time() - start_time
    print("=" * 80)
    print(f"完成: {output_path}")
    print(f"耗时: {elapsed:.1f} 秒")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(description="Precompute lightweight JASPAR Plants PWM score thresholds.")
    parser.add_argument("--samples", type=int, default=100000, help="每个 motif 的随机背景窗口数量")
    parser.add_argument(
        "--background",
        choices=["all", "uniform", "cs_promoter", "fielder_promoter"],
        default="all",
        help="要构建的背景模型",
    )
    parser.add_argument("--seed", type=int, default=123, help="基础随机种子")
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    jaspar_dir = project_root / "data" / "motif_db" / "jaspar_plants"
    pwm_path = jaspar_dir / "jaspar_plants_pwm.json"

    with pwm_path.open("r", encoding="utf-8") as handle:
        motifs = json.load(handle)

    print(f"已读取 motif 数量: {len(motifs)}")

    jobs = []
    if args.background in ("all", "uniform"):
        jobs.append(
            (
                "uniform",
                {base: 0.25 for base in BASES},
                jaspar_dir / "jaspar_background_thresholds_uniform.json",
            )
        )

    cs_background_path = jaspar_dir / "background_cs_promoter_2000.json"
    if args.background in ("all", "cs_promoter"):
        background = load_background_json(cs_background_path)
        if background is None:
            print(f"[INFO] 未找到 CS promoter background，跳过: {cs_background_path}")
        else:
            jobs.append(
                (
                    "cs_promoter",
                    background,
                    jaspar_dir / "jaspar_background_thresholds_cs_promoter.json",
                )
            )

    fielder_background_path = jaspar_dir / "background_fielder_promoter_2000.json"
    if args.background in ("all", "fielder_promoter"):
        background = load_background_json(fielder_background_path)
        if background is None:
            print(f"[INFO] 未找到 Fielder promoter background，跳过: {fielder_background_path}")
        else:
            jobs.append(
                (
                    "fielder_promoter",
                    background,
                    jaspar_dir / "jaspar_background_thresholds_fielder_promoter.json",
                )
            )

    if not jobs:
        raise SystemExit("没有可构建的背景模型。请先运行 scripts/build_promoter_background.py 或选择 --background uniform。")

    total_start = time.time()
    for background_name, background, output_path in jobs:
        build_one_background(
            motifs=motifs,
            background_name=background_name,
            background=background,
            output_path=output_path,
            n_samples=args.samples,
            base_seed=args.seed,
            p_levels=DEFAULT_P_LEVELS,
        )

    print(f"全部完成，总耗时: {time.time() - total_start:.1f} 秒")


if __name__ == "__main__":
    main()
