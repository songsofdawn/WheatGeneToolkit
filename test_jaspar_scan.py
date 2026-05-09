from pathlib import Path

from utils.jaspar_pwm_scan import (
    load_jaspar_pwm,
    load_precomputed_thresholds,
    parse_fasta_or_plain,
    scan_sequences_with_jaspar_thresholds,
)


TEST_FASTA = """>test_promoter
ATGCGTACGTACGTACGTGACGTAGCTAGCTAGCTACGATCGATCGATCGATCGATCACGTGTTGACCTTAA
"""


def main():
    project_root = Path(__file__).resolve().parent
    jaspar_dir = project_root / "data" / "motif_db" / "jaspar_plants"
    pwm_path = jaspar_dir / "jaspar_plants_pwm.json"
    threshold_path = jaspar_dir / "jaspar_background_thresholds_uniform.json"
    output_path = project_root / "test_jaspar_scan_results.csv"

    motifs = load_jaspar_pwm(pwm_path)
    thresholds = load_precomputed_thresholds(threshold_path)
    records = parse_fasta_or_plain(TEST_FASTA)

    print(f"Loaded motifs: {len(motifs)}")
    print(f"Parsed records: {list(records.keys())}")

    if thresholds:
        print(f"Loaded threshold table: {threshold_path}")
        print(f"Threshold samples: {thresholds.get('n_samples')}")
        test_motifs = motifs
    else:
        print(f"Threshold table not found: {threshold_path}")
        print("Run: python scripts/build_jaspar_background.py")
        print("Fallback test uses first 50 motifs without threshold significance.")
        test_motifs = motifs[:50]

    result_df = scan_sequences_with_jaspar_thresholds(
        records=records,
        motifs=test_motifs,
        relative_cutoff=0.90,
        p_level_cutoff=0.01,
        scan_reverse=True,
        max_total_hits=20000,
        precomputed_thresholds=thresholds,
        top_n_per_motif_sequence=3,
    )

    required_columns = {"p_level", "significance_rank", "confidence_level"}
    missing_columns = required_columns - set(result_df.columns)
    if missing_columns:
        raise RuntimeError(f"Missing expected columns: {sorted(missing_columns)}")

    result_df.to_csv(output_path, index=False, encoding="utf-8-sig")
    print(f"Candidate hits: {len(result_df)}")
    if not result_df.empty:
        print(result_df[["matrix_id", "relative_score", "p_level", "confidence_level"]].head().to_string(index=False))
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
