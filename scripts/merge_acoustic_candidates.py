#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_csv(path: Path):
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def as_float(row, key, default=0.0):
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge acoustic-Gorkov candidates from sweep folders."
    )
    parser.add_argument("--root", default="results/week4_final")
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    root = Path(args.root)
    candidate_rows = []
    excitation_rows = []

    for case_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        case_id = case_dir.name
        for row in read_csv(case_dir / "gorkov_candidates.csv"):
            row["acoustic_case_id"] = case_id
            candidate_rows.append(row)
        for row in read_csv(case_dir / "bubble_excitations.csv"):
            row["acoustic_case_id"] = case_id
            excitation_rows.append(row)

    candidate_fieldnames = [
        "acoustic_case_id",
        "candidate_id",
        "node_index",
        "r_m",
        "z_m",
        "pressure_re_pa",
        "pressure_im_pa",
        "pressure_abs_pa",
        "pressure_phase_rad",
        "drive_amplitude_pa",
        "drive_phase_rad",
        "gorkov_potential_j",
        "gorkov_force_abs_n",
        "velocity_abs_m_s",
        "local_potential_minimum",
        "score",
    ]
    excitation_fieldnames = [
        "acoustic_case_id",
        "case_id",
        "frequency_hz",
        "omega_rad_s",
        "bubble_r_m",
        "bubble_z_m",
        "static_pressure_pa",
        "liquid_temperature_k",
        "equilibrium_radius_m",
        "pb_real_pa",
        "pb_imag_pa",
        "pb_abs_pa",
        "pb_argument_rad",
        "drive_amplitude_pa",
        "drive_phase_rad",
    ]

    candidate_rows.sort(key=lambda row: as_float(row, "score"), reverse=True)
    write_csv(root / "all_candidates.csv", candidate_rows, candidate_fieldnames)
    write_csv(root / "bubble_excitations_all.csv", excitation_rows, excitation_fieldnames)
    write_csv(root / "best_acoustic_candidates.csv", candidate_rows[: args.top], candidate_fieldnames)

    print(f"Merged candidates: {len(candidate_rows)}")
    print(f"Merged bubble excitations: {len(excitation_rows)}")
    print(f"Top candidates written to {root / 'best_acoustic_candidates.csv'}")


if __name__ == "__main__":
    main()
