#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_csv(path: Path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows, fieldnames):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Join acoustic candidate table with bubble solver metrics by case_id."
    )
    parser.add_argument("--acoustic", default="results/week4_final/bubble_excitations_all.csv")
    parser.add_argument("--bubble", required=True)
    parser.add_argument("--output", default="results/week4_final/joined_acoustic_bubble_results.csv")
    args = parser.parse_args()

    acoustic_rows = read_csv(Path(args.acoustic))
    bubble_rows = read_csv(Path(args.bubble))
    bubble_by_case = {row.get("case_id", ""): row for row in bubble_rows}

    joined = []
    fieldnames = []
    for row in acoustic_rows:
        case_id = row.get("case_id", "")
        merged = dict(row)
        merged.update(bubble_by_case.get(case_id, {}))
        joined.append(merged)
        for key in merged:
            if key not in fieldnames:
                fieldnames.append(key)

    write_csv(Path(args.output), joined, fieldnames)
    print(f"Joined rows: {len(joined)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
