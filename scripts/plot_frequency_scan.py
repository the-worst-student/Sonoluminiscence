#!/usr/bin/env python3
import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def ReadRows(PathValue):
    with open(PathValue, newline="") as File:
        return list(csv.DictReader(File))


def ToFloat(Value):
    if Value is None or Value == "":
        return math.nan
    try:
        return float(Value)
    except ValueError:
        return math.nan


def SortByFrequency(Rows):
    return sorted(Rows, key=lambda Row: ToFloat(Row.get("frequency_hz")))


def SaveLinePlot(OutputPath, X, Series, Title, XLabel, YLabel):
    plt.figure(figsize=(10, 6))
    for Label, Y in Series:
        plt.plot(X, Y, marker="o", label=Label)
    plt.title(Title)
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OutputPath, dpi=180)
    plt.close()


def PrintTopByColumn(Rows, ColumnName, Title, Count=10):
    Sorted = sorted(Rows, key=lambda Row: ToFloat(Row.get(ColumnName)), reverse=True)
    print("\n" + Title)
    print("case_id,frequency_hz," + ColumnName)
    for Row in Sorted[:Count]:
        CaseId = Row.get("case_id") or Row.get("best_candidate_id") or ""
        print(f"{CaseId},{Row.get('frequency_hz', '')},{Row.get(ColumnName, '')}")


def PlotFrequencyScan(OutputDir):
    SummaryPath = OutputDir / "frequency_scan_summary.csv"
    if not SummaryPath.exists():
        raise FileNotFoundError(f"No frequency_scan_summary.csv in {OutputDir}")

    Rows = SortByFrequency(ReadRows(SummaryPath))
    GraphDir = OutputDir / "plots"
    GraphDir.mkdir(parents=True, exist_ok=True)

    Frequencies = [ToFloat(Row["frequency_hz"]) for Row in Rows]

    SaveLinePlot(
        GraphDir / "01_pressure_vs_frequency.png",
        Frequencies,
        [
            ("manual bubble point |Pb|, Pa", [ToFloat(Row.get("manual_pb_abs_pa")) for Row in Rows]),
            ("best Gorkov candidate |Pb|, Pa", [ToFloat(Row.get("best_pressure_abs_pa")) for Row in Rows]),
            ("max field |P|, Pa", [ToFloat(Row.get("max_pressure_abs_pa")) for Row in Rows]),
        ],
        "Acoustic pressure response over frequency scan",
        "frequency, Hz",
        "pressure amplitude, Pa",
    )

    SaveLinePlot(
        GraphDir / "02_candidate_score_vs_frequency.png",
        Frequencies,
        [("best candidate score", [ToFloat(Row.get("best_score")) for Row in Rows])],
        "Best candidate score over frequency scan",
        "frequency, Hz",
        "score",
    )

    SaveLinePlot(
        GraphDir / "03_gorkov_curvature_vs_frequency.png",
        Frequencies,
        [
            ("lambda_min Hessian(U_G), J/m^2", [ToFloat(Row.get("best_gorkov_curvature_lambda_min_j_m2")) for Row in Rows]),
            ("lambda_max Hessian(U_G), J/m^2", [ToFloat(Row.get("best_gorkov_curvature_lambda_max_j_m2")) for Row in Rows]),
            ("trace Hessian(U_G), J/m^2", [ToFloat(Row.get("best_gorkov_curvature_trace_j_m2")) for Row in Rows]),
        ],
        "Gorkov potential curvature at best candidate",
        "frequency, Hz",
        "curvature, J/m^2",
    )

    SaveLinePlot(
        GraphDir / "04_gorkov_force_vs_frequency.png",
        Frequencies,
        [("|F_G|, N", [ToFloat(Row.get("best_gorkov_force_abs_n")) for Row in Rows])],
        "Gorkov force magnitude at best candidate",
        "frequency, Hz",
        "force, N",
    )

    PrintTopByColumn(Rows, "best_pressure_abs_pa", "Top acoustic candidates by pressure amplitude")
    PrintTopByColumn(Rows, "best_score", "Top acoustic candidates by combined score")
    PrintTopByColumn(Rows, "best_gorkov_curvature_lambda_min_j_m2", "Top acoustic candidates by minimum positive curvature")

    print(f"\nPlots written to: {GraphDir}")


def PlotBubbleResults(OutputDir, BubbleResultsPath):
    if BubbleResultsPath is None or not BubbleResultsPath.exists():
        return

    Rows = SortByFrequency(ReadRows(BubbleResultsPath))
    if not Rows:
        return

    GraphDir = OutputDir / "plots"
    GraphDir.mkdir(parents=True, exist_ok=True)
    Frequencies = [ToFloat(Row.get("frequency_hz")) for Row in Rows]

    SaveLinePlot(
        GraphDir / "05_bubble_temperature_vs_frequency.png",
        Frequencies,
        [("T_max, K", [ToFloat(Row.get("T_max_k")) for Row in Rows])],
        "Bubble maximum temperature over candidate runs",
        "frequency, Hz",
        "T_max, K",
    )

    SaveLinePlot(
        GraphDir / "06_bubble_collapse_metrics_vs_frequency.png",
        Frequencies,
        [
            ("K_R = R0/R_min", [ToFloat(Row.get("K_R")) for Row in Rows]),
            ("K_exp = R_max/R0", [ToFloat(Row.get("K_exp")) for Row in Rows]),
            ("K_range = R_max/R_min", [ToFloat(Row.get("K_range")) for Row in Rows]),
        ],
        "Bubble collapse ratios over candidate runs",
        "frequency, Hz",
        "dimensionless ratio",
    )

    SaveLinePlot(
        GraphDir / "07_bubble_pressure_mach_vs_frequency.png",
        Frequencies,
        [
            ("p_g_max / 1e8", [ToFloat(Row.get("p_g_max_pa")) / 1.0e8 for Row in Rows]),
            ("M_l", [ToFloat(Row.get("M_l")) for Row in Rows]),
        ],
        "Gas pressure and liquid Mach diagnostic over candidate runs",
        "frequency, Hz",
        "scaled value",
    )

    PrintTopByColumn(Rows, "T_max_k", "Top bubble runs by T_max")
    PrintTopByColumn(Rows, "K_R", "Top bubble runs by compression ratio")


def Main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/plot_frequency_scan.py <scan_output_dir> [bubble_results_all.csv]")
        return 1

    OutputDir = Path(sys.argv[1])
    BubbleResultsPath = Path(sys.argv[2]) if len(sys.argv) >= 3 else OutputDir / "bubble_results_all.csv"

    PlotFrequencyScan(OutputDir)
    PlotBubbleResults(OutputDir, BubbleResultsPath)
    return 0


if __name__ == "__main__":
    raise SystemExit(Main())
