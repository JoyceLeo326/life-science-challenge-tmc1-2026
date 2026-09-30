"""Build current library and five-split historical figures from frozen source tables."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("PACKAGE_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "figures"
LIB = ROOT / "plot_inputs" / "CURRENT_MANIFEST.json"
SPLITS = ROOT / "plot_inputs" / "benchmark_five_splits.csv"
AUDIT_LIB = ROOT / "05_independent_audit" / "final_library_independent_audit.json"
AUDIT_BENCH = ROOT / "05_independent_audit" / "benchmark_multisplit_independent_audit.json"
PREFLIGHT = ROOT / "plot_inputs" / "ranking_formal_v2_16698" / "preflight_v2.json"
PREDICTIONS = ROOT / "plot_inputs" / "ranking_formal_v2_16698" / "all_predictions.csv"
AUDIT_DESIGN = ROOT / "05_independent_audit" / "design40_job_independent_qc_summary.json"
AUDIT_DESIGN_PROV = ROOT / "05_independent_audit" / "design40_provenance_independent_audit.json"
DESIGN_CONTACTS = ROOT / "05_independent_audit" / "design40_pose_region_contacts.csv"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(fig, name: str, caption: str, data_rows: list[dict], sources: list[Path]) -> None:
    prefix = OUT / name
    with prefix.with_name(name + "_data.csv").open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(data_rows[0]))
        writer.writeheader()
        writer.writerows(data_rows)
    for suffix in (".svg", ".pdf", ".png"):
        fig.savefig(prefix.with_suffix(suffix), dpi=350, facecolor="white")
    plt.close(fig)
    prefix.with_name(name + "_caption.md").write_text(caption, encoding="utf-8")
    prefix.with_name(name + "_sources.json").write_text(
        json.dumps({"sources": [{"path": path.relative_to(ROOT).as_posix(), "sha256": sha(path)} for path in sources],
                    "script": Path(__file__).resolve().relative_to(ROOT).as_posix(), "script_sha256": sha(Path(__file__)),
                    "visual_checked": False}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def library_flow() -> None:
    current = json.loads(LIB.read_text(encoding="utf-8"))
    preflight = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    with PREDICTIONS.open("r", encoding="utf-8-sig", newline="") as fh:
        predictions = list(csv.DictReader(fh))
    eligible_count = sum(row["docking_eligibility"] == "eligible" for row in predictions)
    counts = current["counts"]
    stages = [
        ("Source records", counts["source_records"], "source records"),
        ("Distinct source IDs", counts["unique_source_ids"], "source IDs"),
        ("Deduplicated parents", counts["standardized_unique_parents"], "parent structures"),
        ("Property eligible", counts["property_filtered_unique_parents"], "parent structures"),
        ("Forward unlabeled", counts["forward_unlabeled_unique_parents"], "parent structures"),
        ("Docking eligible", eligible_count, "parent structures after chemistry gate"),
    ]
    if len(predictions) != 16698 or eligible_count != 9098 or preflight["physical_unique"] != 278:
        raise ValueError("formal prediction/freeze count changed; re-audit before plotting")
    rows = [{"stage_order": i + 1, "stage": name, "count": count, "count_unit": unit}
            for i, (name, count, unit) in enumerate(stages)]
    fig, ax = plt.subplots(figsize=(10, 5.3), layout="constrained")
    colors = ["#9AB7C9", "#6F9AB5", "#3F7EA4", "#226C98", "#D47A42", "#A45978"]
    bars = ax.barh([x[0] for x in stages][::-1], [x[1] for x in stages][::-1], color=colors[::-1])
    for bar, (_, count, _) in zip(bars, stages[::-1]):
        ax.text(count + 280, bar.get_y() + bar.get_height()/2, f"{count:,}", va="center", fontsize=11)
    ax.set_xlim(0, 27000)
    ax.set_xlabel("Count (stage-specific unit shown in caption)", fontsize=11)
    ax.set_title("Public library: frozen source-to-forward flow", fontsize=14, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    caption = ("# F1 正式公共库逐级分母\n\n"
               "ChEMBL来源记录23,475条（其中第一轮只读3,475、本轮API新增20,000）；来源ID去重23,403；"
               "标准化及母体去重21,356个；性质合格17,645个；排除历史来源身份后前瞻未标注16,698个；"
               "其中经03正式v2化学/立体门禁可作业9,098个。"
               "前两级单位为记录或来源ID，后四级为去盐/互变/立体合并后的独立母体。"
               "来源为02不可变快照 `CURRENT_MANIFEST.json` 与03 `preflight_v2.json`，05已独立复核。"
               "安全二维SDF 17,632条仅是可无损表示的交付格式记录数，不是性质合格母体分母。"
               "本图不表示对接成功或生物活性。\n")
    save(fig, "F1_正式公共库分母", caption, rows, [LIB, PREFLIGHT, PREDICTIONS, AUDIT_LIB])


def five_splits() -> None:
    with SPLITS.open("r", encoding="utf-8-sig", newline="") as fh:
        source = list(csv.DictReader(fh))
    methods = [("Random", "random_hits_at96"), ("MW only", "MW_only_hits_at96"),
               ("Ridge", "ridge_hits_at96"), ("ExtraTrees static", "AI_static_hits_at96"),
               ("ExtraTrees active", "AI_active_hits_at96")]
    rows = []
    matrix = []
    for label, key in methods:
        recalls = []
        for record in source:
            denominator = int(record["heldout_top10pct_n"])
            hits = int(record[key])
            recalls.append(hits / denominator)
            rows.append({"split_seed": record["seed"], "train_n": record["train_n"],
                         "heldout_n": record["heldout_n"], "top10pct_n": denominator,
                         "budget_labels": 96, "method": label, "hits": hits,
                         "recall": round(hits / denominator, 9)})
        matrix.append(recalls)
    fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
    heat = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap="YlGnBu")
    ax.set_xticks(range(len(source)), [f"Split {i}\nn={r['heldout_top10pct_n']}" for i, r in enumerate(source, 1)])
    ax.set_yticks(range(len(methods)), [label for label, _ in methods])
    for y, (label, key) in enumerate(methods):
        for x, record in enumerate(source):
            hits = int(record[key]); total = int(record["heldout_top10pct_n"])
            ax.text(x, y, f"{hits}/{total}", ha="center", va="center",
                    color="white" if hits / total >= 0.7 else "#23313B", fontsize=10, weight="bold")
    bar = fig.colorbar(heat, ax=ax, shrink=0.85, pad=0.02)
    bar.set_label("Retrospective top-decile recall")
    ax.set_xlabel("Scaffold-group heldout split; n = top-decile labels", fontsize=11)
    ax.set_title("Historical Vina-score surrogate: 96-label budget per split", fontsize=14, weight="bold")
    ax.tick_params(labelsize=10)
    caption = ("# F2 五次骨架留出固定预算历史回放\n\n"
               "935条首轮质控Vina e8分数作为代理标签；每个Murcko骨架留出切分以96次历史标签查询为预算。"
               "纵轴为该切分已查询Top10%分数标签数/该切分Top10%总数，五次分母分别为36、25、27、21、20。"
               "各切分独立分母不可把命中数直接加总称为同一总体召回。静态ExtraTrees对Ridge胜/平/负0/3/2；"
               "主动ExtraTrees为1/2/2。Top10%由各留出真实对接分数事后确定；这不是前瞻实验或活性预测。"
               "源表、逐次人数及耗时在 `F2_五切分固定预算_data.csv` 和03原始 `benchmark_five_splits.csv`。\n")
    save(fig, "F2_五切分固定预算", caption, rows, [SPLITS, AUDIT_BENCH])


def design_funnel() -> None:
    qc = json.loads(AUDIT_DESIGN.read_text(encoding="utf-8"))
    stages = [
        ("Fixed for docking", 40, "fixed designs"),
        ("Docking completed", int(qc["jobs_in_ledger"]), "actual jobs"),
        ("Technical QC passed", int(qc["independent_qc_pass"]), "jobs"),
        ("Middle contact ≤4 Å", int(qc["middle_contact_4A_among_qc_pass"]), "QC-passed jobs"),
        ("QC + score + Middle", 0, "joint-criterion jobs"),
    ]
    if [x[1] for x in stages[:4]] != [40, 30, 29, 14]:
        raise ValueError("design evidence changed; review before plotting")
    rows = [{"stage_order": i+1, "stage": label, "count": count, "denominator_kind": kind}
            for i, (label, count, kind) in enumerate(stages)]
    fig, ax = plt.subplots(figsize=(10, 5.1), layout="constrained")
    colors = ["#8AA5BA", "#5E91AE", "#2E769F", "#D09355", "#A94F50"]
    bars = ax.barh([x[0] for x in stages][::-1], [x[1] for x in stages][::-1], color=colors[::-1])
    for bar, (_, count, _) in zip(bars, stages[::-1]):
        ax.text(count + 0.5, bar.get_y() + bar.get_height()/2, str(count), va="center", fontsize=11)
    ax.set_xlim(0, 47)
    ax.set_xlabel("Count (stage-specific denominator)", fontsize=11)
    ax.set_title("Generated designs: technical and site-specific attrition", fontsize=14, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    caption = ("# F5 正式SELFIES生成候选的实际复算与联合门槛\n\n"
               "04固定40个设计物送03统一Vina e8流程；10个制备前阻断，30个产生实际对接作业，"
               "05独立技术QC通过29/30。通过QC的29条中，14条的第一姿势与预定15个中区残基有≤4 Å接触。"
               "冻结分数阈值为旧935标签定义的−8.0796 kcal/mol；30条中2条分数达到阈值，"
               "但均无中区接触，其中1条还越盒，故同时满足技术QC、分数阈值和中区接触为0。"
               "最后一栏是联合交集，不是与上一栏同样的单独过滤规则。此图不代表功能活性或生成新颖性。"
               "原始作业与失败见03设计40账本；05独立QC和区域表见来源记录。\n")
    save(fig, "F5_生成设计联合门槛", caption, rows, [AUDIT_DESIGN, AUDIT_DESIGN_PROV, DESIGN_CONTACTS])


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    library_flow()
    five_splits()
    design_funnel()
