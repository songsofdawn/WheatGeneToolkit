from pathlib import Path

from utils.jaspar_pwm_scan import (
    load_jaspar_pwm,
    parse_fasta_or_plain,
    scan_sequences_with_jaspar,
)


TEST_FASTA = """>test_promoter
ATGCGTACGTACGTACGTGACGTAGCTAGCTAGCTACGATCGATCGATCGATCGATCACGTGTTGACCTTAA
"""


def main():
    project_root = Path(__file__).resolve().parent
    pwm_path = project_root / "data" / "motif_db" / "jaspar_plants" / "jaspar_plants_pwm.json"
    output_path = project_root / "test_jaspar_scan_results.csv"

    motifs = load_jaspar_pwm(pwm_path)
    records = parse_fasta_or_plain(TEST_FASTA)

    print(f"Loaded motifs: {len(motifs)}")
    print(f"Parsed records: {list(records.keys())}")

    result_df = scan_sequences_with_jaspar(
        records=records,
        motifs=motifs,
        cutoff=0.85,
        scan_reverse=True,
        max_total_hits=20000,
    )

    result_df.to_csv(output_path, index=False, encoding="utf-8-sig")
    print(f"Hits: {len(result_df)}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
