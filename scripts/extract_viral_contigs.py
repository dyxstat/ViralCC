#!/usr/bin/env python3
"""
Extract viral (non-provirus) contig names from geNomad TSV.

Input : genomad_output/final_contigs_summary/final_contigs_virus_summary.tsv
Output: viral_contig.txt   (one contig name per line, no header)

Rules:
- Keep rows where topology is NOT a provirus.
- We treat 'Provirus' case-insensitively; anything containing 'provirus' is excluded.
- Also exclude if the seq_name itself contains '|provirus_' (extra safety).

Usage:
  python extract_viral_contigs.py \
      --input genomad_output/final_contigs_summary/final_contigs_virus_summary.tsv \
      --output viral_contig.txt
"""

import argparse
import csv
import sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", "-i",
                    default="genomad_output/final_contigs_summary/final_contigs_virus_summary.tsv",
                    help="Path to geNomad virus summary TSV")
    ap.add_argument("--output", "-o", default="viral_contig.txt",
                    help="Path to write contig names (one per line)")
    args = ap.parse_args()

    try:
        with open(args.input, "r", encoding="utf-8", errors="replace", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            required = {"seq_name", "topology"}
            if not required.issubset(reader.fieldnames or []):
                sys.exit(f"ERROR: TSV missing columns {required - set(reader.fieldnames or [])}. "
                         f"Got columns: {reader.fieldnames}")

            kept = []
            for row in reader:
                name = (row.get("seq_name") or "").strip()
                topo = (row.get("topology") or "").strip().lower()

                if not name:
                    continue
                # Exclude proviruses/prophage: topology contains 'provirus'
                if "provirus" in topo:
                    continue
                # Extra guard: some seq_names carry a provirus tag
                if "|provirus_" in name.lower():
                    continue

                kept.append(name)

        with open(args.output, "w", encoding="utf-8") as out:
            if kept:
                out.write("\n".join(kept) + "\n")

        print(f"Done. Wrote {len(kept)} contig names to {args.output}")

    except FileNotFoundError:
        sys.exit(f"ERROR: Input file not found: {args.input}")
    except Exception as e:
        sys.exit(f"ERROR: {e}")

if __name__ == "__main__":
    main()
