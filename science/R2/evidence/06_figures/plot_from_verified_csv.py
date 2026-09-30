"""Render F1-F5 only from nonempty, validated second-round CSV data.

Usage:
  python plot_from_verified_csv.py F1 path/to/stages.csv output/prefix
  python plot_from_verified_csv.py F2 path/to/methods.csv output/prefix
  python plot_from_verified_csv.py F3 path/to/heldout.csv output/prefix
  python plot_from_verified_csv.py F4 path/to/conditions.csv output/prefix
  python plot_from_verified_csv.py F5 path/to/paired.csv output/prefix

F1 columns: stage_order,stage,count,source_path
F2 columns: method,budget_dockings,hit_count,hit_definition,source_path
F3 columns: compound_id,predicted_kcal_mol,actual_kcal_mol,split,scaffold_disjoint,source_path
F4 columns: candidate_id,condition_id,score_kcal_mol,source_path
F5 columns: pair_id,candidate_class,score_kcal_mol,passes_constraints,source_path

F2 methods are compared only within one fixed budget and hit definition. The
caller must verify that the compared methods used the same eligible pool.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone, timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REQUIRED = {
    "F1": {"stage_order", "stage", "count", "source_path"},
    "F2": {"method", "budget_dockings", "hit_count", "hit_definition", "source_path"},
    "F3": {"compound_id", "predicted_kcal_mol", "actual_kcal_mol", "split", "scaffold_disjoint", "source_path"},
    "F4": {"candidate_id", "condition_id", "score_kcal_mol", "source_path"},
    "F5": {"pair_id", "candidate_class", "score_kcal_mol", "passes_constraints", "source_path"},
}


def read_rows(path: Path, figure_id: str) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        missing = REQUIRED[figure_id] - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("source table is empty; no figure rendered")
    for index, row in enumerate(rows, 2):
        if not all(row.get(key, "").strip() for key in REQUIRED[figure_id]):
            raise ValueError(f"row {index} has an empty required field")
        source = Path(row["source_path"])
        if not source.is_file():
            raise ValueError(f"row {index} missing evidence file: {source}")
    return rows


def draw_f1(ax, rows: list[dict[str, str]]) -> str:
    ordered = sorted(rows, key=lambda row: int(row["stage_order"]))
    stages = [row["stage"] for row in ordered]
    counts = [int(row["count"]) for row in ordered]
    if len(stages) != len(set(stages)) or any(n < 0 for n in counts):
        raise ValueError("stage names must be unique and counts nonnegative")
    if any(later > earlier for earlier, later in zip(counts, counts[1:])):
        raise ValueError("flow counts increased; stage definition must be checked")
    bars = ax.barh(stages[::-1], counts[::-1], color="#2679A8")
    for bar, count in zip(bars, counts[::-1]):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2, f"  {count:,}", va="center")
    ax.set_xlim(0, max(counts) * 1.19 if max(counts) else 1)
    ax.set_xlabel("Records / unique source IDs / unique parents (n; see stage)")
    ax.set_title("Public library processing flow")
    return f"Observed counts across {len(rows)} ordered stages; entry denominator: {counts[0]:,} source records. Stage labels distinguish source IDs from deduplicated parents."


def draw_f2(ax, rows: list[dict[str, str]]) -> str:
    budgets = {int(row["budget_dockings"]) for row in rows}
    definitions = {row["hit_definition"] for row in rows}
    if len(budgets) != 1 or len(definitions) != 1:
        raise ValueError("F2 requires one shared docking budget and hit definition")
    budget = budgets.pop()
    if budget <= 0 or len(rows) < 2:
        raise ValueError("F2 requires a positive budget and at least two methods")
    names = [row["method"] for row in rows]
    hits = [int(row["hit_count"]) for row in rows]
    if len(names) != len(set(names)) or any(n < 0 or n > budget for n in hits):
        raise ValueError("duplicate methods or invalid hit counts")
    palette = ["#2679A8", "#D28C32", "#659B73", "#A8789A"]
    bars = ax.barh(names[::-1], hits[::-1], color=[palette[i % len(palette)] for i in range(len(rows))][::-1])
    for bar, hit in zip(bars, hits[::-1]):
        ax.text(hit, bar.get_y() + bar.get_height() / 2, f"  {hit}", va="center")
    ax.set_xlabel("Candidates meeting predefined docking criterion (n)")
    ax.set_xlim(0, max(hits) * 1.2 + 0.5)
    ax.set_title("Fixed-budget selection comparison")
    return f"Methods share a budget of {budget} actual docking evaluations each. Hit definition: {definitions.pop()}. Outcome is a docking proxy, not biological activity."


def draw_f3(ax, rows: list[dict[str, str]]) -> str:
    if any(row["split"].lower() != "heldout" or row["scaffold_disjoint"].lower() != "yes" for row in rows):
        raise ValueError("F3 accepts only documented scaffold-disjoint heldout rows")
    ids = [row["compound_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate heldout compound_id")
    predicted = [float(row["predicted_kcal_mol"]) for row in rows]
    actual = [float(row["actual_kcal_mol"]) for row in rows]
    low, high = min(predicted + actual), max(predicted + actual)
    ax.scatter(actual, predicted, s=22, alpha=0.65, color="#2679A8")
    ax.plot([low, high], [low, high], color="#666666", linestyle="--", linewidth=1)
    ax.set_xlabel("Actual docking score (kcal/mol)")
    ax.set_ylabel("Predicted docking score (kcal/mol)")
    ax.set_title("Scaffold-disjoint heldout prediction")
    mae = sum(abs(a - p) for a, p in zip(actual, predicted)) / len(rows)
    return f"{len(rows)} unique scaffold-disjoint heldout molecules; mean absolute error {mae:.3f} kcal/mol. This is a docking-score surrogate, not biological activity prediction."


def draw_f4(ax, rows: list[dict[str, str]]) -> str:
    by_candidate: dict[str, list[tuple[str, float]]] = defaultdict(list)
    seen = set()
    for row in rows:
        key = row["candidate_id"], row["condition_id"]
        if key in seen:
            raise ValueError(f"duplicate candidate-condition: {key}")
        seen.add(key)
        by_candidate[key[0]].append((key[1], float(row["score_kcal_mol"])))
    if not by_candidate or any(len(values) < 2 for values in by_candidate.values()):
        raise ValueError("F4 requires at least two observed conditions per candidate")
    names = sorted(by_candidate)
    for x, name in enumerate(names):
        scores = [score for _, score in by_candidate[name]]
        ax.scatter([x] * len(scores), scores, color="#2679A8", alpha=0.7, s=25)
        ax.plot([x - 0.18, x + 0.18], [sum(scores) / len(scores)] * 2, color="#D28C32", linewidth=2)
    ax.set_xticks(range(len(names)), names, rotation=45, ha="right")
    ax.set_ylabel("Docking score (kcal/mol)")
    ax.set_title("Scores across explicitly recorded conditions")
    return f"{len(rows)} observed candidate-condition scores from {len(names)} candidates; orange segments show per-candidate arithmetic means. Scores do not establish functional direction."


def draw_f5(ax, rows: list[dict[str, str]]) -> str:
    pairs: dict[str, dict[str, float]] = defaultdict(dict)
    for row in rows:
        kind = row["candidate_class"].lower()
        if kind not in {"parent", "design"} or row["passes_constraints"].lower() != "yes":
            raise ValueError("F5 requires parent/design rows with verified constraints")
        if kind in pairs[row["pair_id"]]:
            raise ValueError(f"duplicate pair/class: {row['pair_id']}/{kind}")
        pairs[row["pair_id"]][kind] = float(row["score_kcal_mol"])
    if not pairs or any(set(values) != {"parent", "design"} for values in pairs.values()):
        raise ValueError("F5 requires exactly one parent and one design per pair")
    names = sorted(pairs)
    for index, name in enumerate(names):
        parent, design = pairs[name]["parent"], pairs[name]["design"]
        ax.plot([index - 0.17, index + 0.17], [parent, design], color="#9AA3AC", linewidth=1.2)
        ax.scatter([index - 0.17], [parent], color="#2679A8", s=30)
        ax.scatter([index + 0.17], [design], color="#D28C32", s=30)
    ax.set_xticks(range(len(names)), names, rotation=45, ha="right")
    ax.set_ylabel("Docking score (kcal/mol)")
    ax.set_title("Same-protocol parent/design comparison")
    return f"{len(names)} parent/design pairs with verified constraints and actual docking scores. Blue: parent; orange: design. A score difference does not establish experimental improvement."


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figure_id", choices=sorted(REQUIRED))
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_prefix", type=Path)
    args = parser.parse_args()
    rows = read_rows(args.input_csv, args.figure_id)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.5, 5.2), layout="constrained")
    caption = {"F1": draw_f1, "F2": draw_f2, "F3": draw_f3, "F4": draw_f4, "F5": draw_f5}[args.figure_id](ax, rows)
    for suffix in (".svg", ".pdf", ".png"):
        fig.savefig(args.output_prefix.with_suffix(suffix), dpi=350, facecolor="white")
    plt.close(fig)
    args.output_prefix.with_name(args.output_prefix.name + "_caption.md").write_text(
        f"# {args.figure_id} 图注\n\n{caption}\n\n源表：`{args.input_csv.resolve()}`。\n",
        encoding="utf-8",
    )
    evidence_paths = sorted({str(Path(row["source_path"]).resolve()) for row in rows})
    sources = {
        "figure_id": args.figure_id,
        "generated_at_cst": datetime.now(timezone(timedelta(hours=8))).isoformat(),
        "data_path": str(args.input_csv.resolve()),
        "data_sha256": hashlib.sha256(args.input_csv.read_bytes()).hexdigest(),
        "evidence_paths": evidence_paths,
        "visual_checked": False,
    }
    args.output_prefix.with_name(args.output_prefix.name + "_sources.json").write_text(
        json.dumps(sources, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
