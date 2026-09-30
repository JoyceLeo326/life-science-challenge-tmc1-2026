"""Redraw F3/F4 only from the verified, already-computed source CSV tables."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

SOURCE_SHA = {
    "actual_parent_product_comparison.csv": "758f4c230f4e75c28788b96a4356d273f13c09d80baa5149d291305b26ea96be",
    "generator_pilot_same_budget_comparison.csv": "c02b8f876157822c39d8c10d7e24b03823d33c2013c45dba5f22bb3cc4ef27a0",
    "generator_pilot_identity_safe_funnel.csv": "a8cdd739ddbba0b952ae530d6d9a728fd725a0ae20614d2683c90a44273c9873",
}


def rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def f3(data: Path, out: Path) -> None:
    results = rows(data)
    if len(results) != 7 or not all(r["same_protocol"] == "True" and r["independent_graph_pose_score_hash_qc"] == "True" for r in results):
        raise ValueError("F3 source must have seven same-protocol, independently QC-passed pairs")
    fig, ax = plt.subplots(figsize=(8.5, 5.6), dpi=160)
    colors = {"shared_middle_contact": "#176A7A", "region_shift_or_no_middle_contact": "#C65A25"}
    for i, r in enumerate(results):
        parent = float(r["parent_score_kcal_mol"])
        design = float(r["design_score_kcal_mol"])
        color = colors[r["region_class"]]
        ax.plot([parent, design], [i, i], color=color, lw=2, alpha=.8)
        ax.scatter([parent], [i], marker="o", s=48, facecolor="white", edgecolor=color, linewidth=1.6, zorder=3)
        ax.scatter([design], [i], marker="D", s=44, color=color, zorder=3)
        ax.annotate(f"{float(r['design_minus_parent_kcal_mol']):+.2f}", ((parent+design)/2, i-.16),
                    ha="center", va="bottom", fontsize=8)
    ax.set_yticks(range(len(results)), [f"{r['parent_name'].title()} / {r['transformation'].replace('_', ' ')}" for r in results])
    ax.set_xlabel("AutoDock Vina e8 score (kcal/mol; lower is more favorable)")
    ax.set_title("Matched parent and enumerated analogues: single same-protocol run")
    ax.invert_xaxis()
    ax.invert_yaxis()
    ax.set_ylim(len(results)-.35, -.75)
    ax.grid(axis="x", color="#d0d8da", lw=.7)
    ax.set_axisbelow(True)
    legend = [Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#176A7A", label="parent"),
              Line2D([0], [0], marker="D", color="none", markerfacecolor="#176A7A", markeredgecolor="#176A7A", label="analogue with Middle contact <=4 Å"),
              Line2D([0], [0], marker="D", color="none", markerfacecolor="#C65A25", markeredgecolor="#C65A25", label="analogue without Middle contact <=4 Å")]
    ax.legend(handles=legend, loc="lower right", fontsize=8, frameon=False)
    fig.text(.02, .01, "Numbers are ΔVina = analogue − parent (kcal/mol). All 9 jobs: same receptor, box, e8 and seed; independent QC passed. Orange = region shift.", fontsize=7.5)
    fig.tight_layout(rect=(0, .035, 1, 1))
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"actual_parent_product_comparison.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def f4(data: Path, funnel: Path, out: Path) -> None:
    values = rows(data)
    strict = {r["method"]: r for r in rows(funnel)}
    labels = ["SMILES char GRU", "SELFIES GRU base", "SELFIES GRU focused"]
    if [r["method"] for r in values] != labels or set(strict) != set(labels):
        raise ValueError("F4 method rows changed")
    if any(int(r["requested_samples"]) != 10000 for r in values):
        raise ValueError("F4 sample denominator changed")
    for r in values:
        if int(r["identity_safe_property_and_alert_pass_unique"]) != int(strict[r["method"]]["identity_safe_property_alert_pass"]):
            raise ValueError(f"F4 identity-safe funnel differs: {r['method']}")
    ticks = ["valid", "unique valid", "property + alerts", "identity-safe", "similarity >=0.15", "similarity >=0.25"]
    fields = ["valid_samples", "unique_valid_samples", "property_and_alert_pass_unique", "identity_safe_property_and_alert_pass_unique", "property_pass_and_similarity_ge_0p15", "property_pass_and_similarity_ge_0p25"]
    fig, ax = plt.subplots(figsize=(9.2, 5.4), dpi=160)
    x = np.arange(len(fields))
    width = .25
    colors = ["#4A6172", "#176A7A", "#C65A25"]
    for i, r in enumerate(values):
        counts = [int(r[f]) for f in fields]
        rects = ax.bar(x+(i-1)*width, counts, width, color=colors[i], label=r["method"])
        for rect, n in zip(rects, counts):
            ax.text(rect.get_x()+rect.get_width()/2, max(n, .8)*1.15, str(n), ha="center", va="bottom", fontsize=7, rotation=0)
    ax.set_yscale("log")
    ax.set_ylim(.7, 30000)
    ax.set_xticks(x, ticks)
    ax.set_ylabel("Molecules per 10,000 samples (log scale)")
    ax.set_title("Generative model pilot: same library and sample budget, with identity gate")
    ax.grid(axis="y", which="major", color="#dbe2e5", lw=.65)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.text(.02, .01, "Identity-safe excludes radicals and unspecified stereochemistry. All use 5,067 ChEMBL molecules; 10 base epochs, +3 focused epochs. No docking labels.", fontsize=7.2)
    fig.tight_layout(rect=(0, .035, 1, 1))
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"generator_pilot_same_budget_comparison.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r2", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    base = a.r2 / "04_design/method_results"
    data3 = base / "actual_parent_product_comparison.csv"
    data4 = base / "generator_pilot_same_budget_comparison.csv"
    funnel = base / "generator_pilot_identity_safe_funnel.csv"
    for path in (data3, data4, funnel):
        if not path.is_file():
            raise FileNotFoundError(f"F3/F4 verified source CSV missing: {path}")
        if sha(path) != SOURCE_SHA[path.name]:
            raise ValueError(f"F3/F4 verified source CSV SHA256 mismatch: {path}")
    a.output.mkdir(parents=True, exist_ok=True)
    f3(data3, a.output)
    f4(data4, funnel, a.output)
    print({"source_sha256": {str(x): sha(x) for x in (data3, data4, funnel)}, "output": str(a.output)})


if __name__ == "__main__":
    main()
