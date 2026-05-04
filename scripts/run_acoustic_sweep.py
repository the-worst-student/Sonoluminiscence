#!/usr/bin/env python3
import argparse
import re
import shutil
import subprocess
from pathlib import Path


def replace_frequency(config_text: str, frequency_hz: float) -> str:
    return re.sub(
        r"(^\s*frequency_hz:\s*)[-+0-9.eE]+",
        rf"\g<1>{frequency_hz}",
        config_text,
        count=1,
        flags=re.MULTILINE,
    )


def run(command, cwd: Path) -> None:
    print("+", " ".join(str(x) for x in command))
    subprocess.run(command, cwd=cwd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run acoustic-Gorkov sweep over selected frequencies."
    )
    parser.add_argument("--build-dir", default="build")
    parser.add_argument("--config", default="configs/base.yaml")
    parser.add_argument("--output-root", default="results/week4_final")
    parser.add_argument(
        "--frequencies",
        default="20000,30000,40000,60000,80000,100000",
        help="Comma-separated frequencies in Hz.",
    )
    parser.add_argument("--sectors", type=int, default=48)
    parser.add_argument("--max-candidates", type=int, default=10)
    args = parser.parse_args()

    project_root = Path.cwd()
    build_dir = (project_root / args.build_dir).resolve()
    config_path = (project_root / args.config).resolve()
    output_root = (project_root / args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    config_text = config_path.read_text()
    frequencies = [float(x) for x in args.frequencies.split(",") if x.strip()]

    if not (build_dir / "build_geometry").exists():
        raise RuntimeError(f"build_geometry not found in {build_dir}")
    if not (build_dir / "solve_acoustics").exists():
        raise RuntimeError(f"solve_acoustics not found in {build_dir}")

    for frequency in frequencies:
        case_id = f"f_{int(round(frequency))}"
        case_dir = output_root / case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        case_config = case_dir / "config.yaml"
        case_config.write_text(replace_frequency(config_text, frequency))

        mesh_path = case_dir / "mesh.msh"
        excitations_path = case_dir / "bubble_excitations.csv"
        slice_path = case_dir / "acoustic_slice.vtk"
        volume_path = case_dir / "acoustic_pseudovolume.vtk"
        field_path = case_dir / "gorkov_field.csv"
        candidates_path = case_dir / "gorkov_candidates.csv"

        run(["./build_geometry", str(case_config), str(mesh_path)], build_dir)
        run([
            "./solve_acoustics",
            str(case_config),
            str(mesh_path),
            str(excitations_path),
            str(slice_path),
            str(volume_path),
            str(args.sectors),
            str(field_path),
            str(candidates_path),
            str(args.max_candidates),
        ], build_dir)

    shutil.copy2(config_path, output_root / "base_config_used.yaml")
    print(f"Sweep finished. Results are in {output_root}")


if __name__ == "__main__":
    main()
