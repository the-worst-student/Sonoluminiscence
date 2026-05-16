#!/usr/bin/env python3

import csv
import re
import subprocess
from pathlib import Path


BASE_CONFIG = Path("configs/base.yaml")
RUN_PIPELINE = Path("build/run_pipeline")
OUTPUT_DIR = Path("results/vessel_frequency_sweep_py")

FREQUENCIES_KHZ = range(50, 201, 10)

VESSELS = {
    "cylinder": """\
vessel:
  type: "cylinder"
  radius_m: 0.05
  height_m: 0.10
""",
    "conical": """\
vessel:
  type: "conical"
  radius_m: 0.05
  height_m: 0.10
  bottom_radius_m: 0.055
  top_radius_m: 0.040
  profile_points: 2
""",
    "barrel": """\
vessel:
  type: "barrel"
  radius_m: 0.05
  height_m: 0.10
  bottom_radius_m: 0.045
  top_radius_m: 0.045
  bulge_m: 0.012
  profile_points: 50
""",
    "hourglass": """\
vessel:
  type: "hourglass"
  radius_m: 0.05
  height_m: 0.10
  bottom_radius_m: 0.055
  top_radius_m: 0.055
  neck_m: 0.012
  profile_points: 50
""",
}


def replace_vessel_block(text, new_vessel_block):
    return re.sub(
        r"  vessel:\n(?:    .*\n)+",
        "  " + new_vessel_block.replace("\n", "\n  ").rstrip() + "\n",
        text,
        count=1,
        )


def replace_frequency(text, frequency_hz):
    return re.sub(
        r"frequency_hz:\s*[-+0-9.eE]+",
        f"frequency_hz: {frequency_hz}",
        text,
        count=1,
    )


def get_value(text, pattern):
    match = re.search(pattern, text)
    if match:
        return match.group(1)

    return ""


def parse_output(text):
    return {
        "pb_abs_pa": get_value(text, r"Pb abs \[Pa\]:\s*([-+0-9.eE]+)"),
        "pb_phase_rad": get_value(text, r"Pb phase \[rad\]:\s*([-+0-9.eE]+)"),
        "min_radius_m": get_value(text, r"Min radius \[m\]:\s*([-+0-9.eE]+)"),
        "max_temperature_k": get_value(text, r"Max gas temperature \[K\]:\s*([-+0-9.eE]+)"),
        "temperature_at_collapse_k": get_value(text, r"Temperature at collapse \[K\]:\s*([-+0-9.eE]+)"),
        "compression_ratio": get_value(text, r"Compression ratio R0/Rmin:\s*([-+0-9.eE]+)"),
        "potential_luminescence": get_value(text, r"Potential luminescence:\s*(\w+)"),
    }


def main():
    base_text = BASE_CONFIG.read_text()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []

    for vessel_name, vessel_block in VESSELS.items():
        for frequency_khz in FREQUENCIES_KHZ:
            frequency_hz = frequency_khz * 1000.0

            case_dir = OUTPUT_DIR / vessel_name / f"{frequency_khz:03d}_khz"
            case_dir.mkdir(parents=True, exist_ok=True)

            config_text = replace_vessel_block(base_text, vessel_block)
            config_text = replace_frequency(config_text, frequency_hz)

            config_path = case_dir / "config.yaml"
            log_path = case_dir / "run.log"

            config_path.write_text(config_text)

            result = subprocess.run(
                [str(RUN_PIPELINE.resolve()), str(config_path.resolve())],
                cwd=case_dir,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            log_path.write_text(result.stdout)

            metrics = parse_output(result.stdout)

            row = {
                "vessel": vessel_name,
                "frequency_khz": frequency_khz,
                "frequency_hz": frequency_hz,
                "status": result.returncode,
                **metrics,
                "case_dir": str(case_dir),
            }

            rows.append(row)

            print(
                vessel_name,
                frequency_khz,
                "kHz",
                "status =",
                result.returncode,
                "Pb_abs =",
                metrics["pb_abs_pa"],
                "T_max =",
                metrics["max_temperature_k"],
                "K",
                "T_collapse =",
                metrics["temperature_at_collapse_k"],
                "K",
            )

    csv_path = OUTPUT_DIR / "summary.csv"

    with csv_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    successful_rows = [row for row in rows if row["status"] == 0 and row["pb_abs_pa"]]

    if successful_rows:
        best = max(successful_rows, key=lambda row: float(row["pb_abs_pa"]))

        print()
        print("Best by Pb_abs:")
        print("vessel:", best["vessel"])
        print("frequency_khz:", best["frequency_khz"])
        print("Pb_abs_pa:", best["pb_abs_pa"])
        print("T_max_K:", best["max_temperature_k"])
        print("T_collapse_K:", best["temperature_at_collapse_k"])

    print()
    print("CSV written to:", csv_path)


if __name__ == "__main__":
    main()