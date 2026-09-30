"""Plot the frozen three-arm forward comparison from independently audited tables."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(os.environ.get("PACKAGE_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "figures"
SRC = ROOT / "data"
SUMMARY = SRC / "forward_arm_summary.csv"
STRATA = SRC / "forward_strata.csv"
AUDIT = ROOT / "05_independent_audit" / "forward_evaluation_independent_audit.json"
STEM = "F8_正式前瞻三臂固定预算"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    rows = load_csv(SUMMARY)
    strata = load_csv(STRATA)
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    assert audit["all_checks_pass"] is True
    assert audit["physical_unique_queries"] == 319
    assert {r["arm"] for r in rows} == {"random", "descriptor_ridge", "AI"}
    assert all(int(r["attempted"]) == int(r["qc_pass"]) == 120 for r in rows)
    assert sha(SUMMARY) == "77ff7e7e6b70b80e7309983775b0b45cf15235712fbc95d4fa1989bafafa712b"
    by_arm = {r["arm"]: r for r in rows}
    assert [(int(by_arm[a]["favorable_score"]), int(by_arm[a]["middle_contact_4A"]),
             int(by_arm[a]["favorable_and_middle_contact"]))
            for a in ("random", "descriptor_ridge", "AI")] == [(36, 91, 29), (86, 119, 85), (115, 107, 103)]
    assert [int(r["attempted"]) for r in strata if r["stratum"] == "mw_below_350"] == [39, 0, 0]

    labels = ["Random", "Descriptor Ridge", "Active AI"]
    arms = ["random", "descriptor_ridge", "AI"]
    endpoints = [
        ("Favorable score", "favorable_score", "#35688f"),
        ("Middle contact", "middle_contact_4A", "#c77b46"),
        ("Both", "favorable_and_middle_contact", "#4b8c7a"),
    ]
    x = np.arange(3)
    width = 0.235
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    for j, (label, key, color) in enumerate(endpoints):
        values = [int(by_arm[a][key]) for a in arms]
        bars = ax.bar(x + (j-1)*width, values, width=width*0.93, label=label, color=color)
        for bar, value in zip(bars, values):
            ax.annotate(f"{value}", (bar.get_x() + bar.get_width()/2, value),
                        xytext=(0, 4), textcoords="offset points", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(x, labels, fontsize=11)
    ax.set_ylim(0, 132)
    ax.set_yticks([0, 30, 60, 90, 120], ["0", "25", "50", "75", "100"])
    ax.set_ylabel("Share of 120 frozen attempts (%)", fontsize=11)
    ax.set_title("Prospective docking screen: three frozen 120-attempt arms", fontsize=14, weight="bold", pad=12)
    ax.legend(ncol=3, frameon=False, loc="upper left", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.18)
    ax.set_axisbelow(True)
    fig.subplots_adjust(bottom=0.20, top=0.84, left=0.10, right=0.98)
    fig.text(0.01, 0.015,
             "Counts above bars; 120/120 independent technical QC per arm. 319 unique physical jobs across 360 arm attempts.",
             fontsize=9, color="#444444")
    for ext in ("svg", "pdf", "png"):
        fig.savefig(OUT / f"{STEM}.{ext}", dpi=350, facecolor="white")
    plt.close(fig)

    out_rows = []
    for label, arm in zip(labels, arms):
        r = by_arm[arm]
        out_rows.append({
            "arm": arm, "label": label, "attempted": 120, "independent_qc_pass": 120,
            "favorable_score": int(r["favorable_score"]),
            "middle_contact_4A": int(r["middle_contact_4A"]),
            "favorable_and_middle_contact": int(r["favorable_and_middle_contact"]),
            "mw_below_350_attempts": next(int(s["attempted"]) for s in strata if s["arm"] == arm and s["stratum"] == "mw_below_350"),
        })
    with (OUT / f"{STEM}_data.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)
    (OUT / f"{STEM}_caption.md").write_text(
        "# F8 正式前瞻三臂固定预算比较\n\n"
        "冻结随机、理化描述符Ridge、主动AI各120次候选尝试，05独立技术QC均120/120；"
        "三臂合计360名额对应319个去重物理对接作业。蓝色为Vina e8分数≤历史935标签预先冻结的第10百分位阈值"
        "−8.07960033416748 kcal/mol；橙色为首姿势与预设孔道中区4 Å接触；绿色为两条件同时满足。"
        "柱顶数字为尝试数，纵轴为固定120名额中的百分比；单次正式实验没有重复运行误差条。"
        "图形只比较本次静态对接代理筛选。随机臂有39/120个MW<350，Ridge和AI均0/120，"
        "且骨架、相似度与电荷分布不同，不能把臂间差异解释为纯模型因果效应或实验活性优势。"
        "历史五次骨架留出回放中主动AI没有稳定超过Ridge，与本次前瞻固定预算结果分别报告。"
        "两项电荷参数警告保留主分析并另作视作失败敏感性。\n",
        encoding="utf-8",
    )
    sources = [SUMMARY, STRATA, SRC / "forward_charge_sensitivity.csv", AUDIT]
    (OUT / f"{STEM}_sources.json").write_text(
        json.dumps({"sources": [{"path": p.relative_to(ROOT).as_posix(), "sha256": sha(p)} for p in sources],
                    "script": Path(__file__).resolve().relative_to(ROOT).as_posix(), "script_sha256": sha(Path(__file__)),
                    "visual_checked": False}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
