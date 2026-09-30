"""Register figures after human visual inspection; does not claim scientific approval."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path
import os

ROOT = Path(os.environ.get("PACKAGE_ROOT", Path(__file__).resolve().parents[1]))
OUT = ROOT / "06_交付整合" / "科研图表"
DESIGN = ROOT / "04_分子设计"
STRUCTURE = ROOT / "09_结构图"
NOW = datetime.now(timezone(timedelta(hours=8))).isoformat()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


rows = []
for figure_id, name, title, unit, denominator in [
    ("F1", "F1_正式公共库分母", "正式公共库逐级分母", "记录/来源ID/去重母体", "23,475来源记录"),
    ("F2", "F2_五切分固定预算", "历史五切分固定预算回放", "召回率", "各切分Top10%标签36/25/27/21/20；每次96查询"),
    ("F5", "F5_生成设计联合门槛", "生成设计实际复算联合门槛", "分子/作业数", "40固定生成结构"),
]:
    sources = OUT / f"{name}_sources.json"
    obj = json.loads(sources.read_text(encoding="utf-8"))
    obj["visual_checked"] = True
    obj["visual_checked_at_cst"] = NOW
    obj["visual_note"] = "PNG opened and inspected for readable labels, counts and clipping."
    sources.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    rows.append({"figure_id": figure_id, "title": title, "data_path": str(OUT / f"{name}_data.csv"),
                 "script_path": str(OUT / "build_verified_stage_figures.py"),
                 "caption_path": str(OUT / f"{name}_caption.md"), "sources_path": str(sources),
                 "vector_path": str(OUT / f"{name}.svg"), "png_path": str(OUT / f"{name}.png"),
                 "unit": unit, "denominator": denominator, "visual_checked": "yes",
                 "review_status": "stage_research_figure", "notes": "PNG and SVG/PDF generated from frozen real source; not final forward result"})

for figure_id, stem, title, caption, source_table, script, unit, denom in [
    ("F3", "actual_parent_product_comparison", "规则设计与亲本同协议复算", "F3_规则父子配对_caption.md",
     "actual_parent_product_comparison.csv", "analyze_parent_product.py", "kcal/mol", "7对设计/母体"),
    ("F4", "generator_pilot_same_budget_comparison", "生成模型同采样预算先导", "F4_生成先导同预算_caption.md",
     "generator_pilot_same_budget_comparison.csv", "compare_generator_pilots.py", "每10000次采样计数", "各模型10000次采样；首批5067库"),
]:
    sources = OUT / f"{figure_id}_复用来源.json"
    files = [DESIGN / f"{stem}.{ext}" for ext in ("csv", "svg", "pdf", "png")]
    files += [DESIGN / script]
    sources.write_text(json.dumps({"reused_from": "04_分子设计", "files": [{"path": str(p), "sha256": sha(p)} for p in files],
                                  "visual_checked": True, "visual_checked_at_cst": NOW}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    rows.append({"figure_id": figure_id, "title": title, "data_path": str(DESIGN / source_table),
                 "script_path": str(DESIGN / script), "caption_path": str(OUT / caption),
                 "sources_path": str(sources), "vector_path": str(DESIGN / f"{stem}.svg"),
                 "png_path": str(DESIGN / f"{stem}.png"), "unit": unit, "denominator": denom,
                 "visual_checked": "yes", "review_status": "stage_research_figure",
                 "notes": "Reused 04 figure after PNG visual review; PDF also available"})

structure_sources = OUT / "F6_复用来源.json"
structure_files = [STRUCTURE / "图09_六链结构与双区域.png", STRUCTURE / "图09_六链结构与双区域.svg",
                   STRUCTURE / "structure_render_manifest.json", STRUCTURE / "render_structure_panels.py",
                   STRUCTURE / "compose_structure_figure.py"]
structure_sources.write_text(json.dumps({"reused_from": "09_结构图", "files": [{"path": str(p), "sha256": sha(p)} for p in structure_files],
                                          "visual_checked": True, "visual_checked_at_cst": NOW}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
rows.append({"figure_id": "F6", "title": "真实坐标六链模型与双区域", "data_path": str(STRUCTURE / "structure_render_manifest.json"),
             "script_path": str(STRUCTURE / "compose_structure_figure.py"),
             "caption_path": str(OUT / "F6_真实坐标结构图_caption.md"), "sources_path": str(structure_sources),
             "vector_path": str(STRUCTURE / "图09_六链结构与双区域.svg"),
             "png_path": str(STRUCTURE / "图09_六链结构与双区域.png"),
             "unit": "Å", "denominator": "1个冻结小鼠六链坐标模型、2个明确区分的搜索盒", "visual_checked": "yes",
             "review_status": "structural_context_only", "notes": "Real PDB coordinates; SVG has embedded rasterized PyMOL 3D panels"})

pose_sources = OUT / "F7_复用来源.json"
pose_files = [STRUCTURE / "图09_亲本设计物区域迁移.png", STRUCTURE / "图09_亲本设计物区域迁移.svg",
              STRUCTURE / "图09_区域迁移数据.csv", STRUCTURE / "pose_shift_render_manifest.json",
              STRUCTURE / "render_pose_shift.py", STRUCTURE / "compose_pose_shift_figure.py"]
pose_sources.write_text(json.dumps({"reused_from": "09_结构图", "files": [{"path": str(p), "sha256": sha(p)} for p in pose_files],
                              "visual_checked": True, "visual_checked_at_cst": NOW,
                              "visual_note": "PNG opened and inspected after independent six-job QC confirmation."},
                             ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
rows.append({"figure_id": "F7", "title": "Trilaciclib去甲基设计物跨种子位姿迁移",
             "data_path": str(STRUCTURE / "图09_区域迁移数据.csv"),
             "script_path": str(STRUCTURE / "compose_pose_shift_figure.py"),
             "caption_path": str(OUT / "F7_亲本设计物位姿迁移_caption.md"),
             "sources_path": str(pose_sources),
             "vector_path": str(STRUCTURE / "图09_亲本设计物区域迁移.svg"),
             "png_path": str(STRUCTURE / "图09_亲本设计物区域迁移.png"),
             "unit": "kcal/mol；Å", "denominator": "3个实际种子、6条作业，05独立QC 6/6",
             "visual_checked": "yes", "review_status": "computational_pose_comparison_only",
             "notes": "Real first-pose coordinates and independent QC; score gain with region shift does not show affinity gain"})

forward_sources = OUT / "F8_正式前瞻三臂固定预算_sources.json"
forward_files = [ROOT / "03_算法与计算" / "formal_forward_v2_final_evaluation" / n
                 for n in ("forward_arm_summary.csv", "forward_strata.csv", "forward_charge_sensitivity.csv")]
forward_files.append(ROOT / "05_独立方法复核" / "forward_evaluation_independent_audit.json")
forward_sources.write_text(json.dumps({"files": [{"path": str(p), "sha256": sha(p)} for p in forward_files],
                                       "script": str(OUT / "plot_forward_three_arms.py"),
                                       "script_sha256": sha(OUT / "plot_forward_three_arms.py"),
                                       "visual_checked": False}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
rows.append({"figure_id": "F8", "title": "正式前瞻三臂固定预算比较",
             "data_path": str(OUT / "F8_正式前瞻三臂固定预算_data.csv"),
             "script_path": str(OUT / "plot_forward_three_arms.py"),
             "caption_path": str(OUT / "F8_正式前瞻三臂固定预算_caption.md"),
             "sources_path": str(forward_sources),
             "vector_path": str(OUT / "F8_正式前瞻三臂固定预算.svg"),
             "png_path": str(OUT / "F8_正式前瞻三臂固定预算.png"),
             "unit": "尝试数/120；%", "denominator": "随机、Ridge、AI各120固定尝试；319个去重物理作业",
             "visual_checked": "no", "review_status": "pending_visual_review",
             "notes": "Independent 05 audit required; static docking proxy; strata/OOD distributions differ"})

with (OUT / "图表登记.csv").open("w", encoding="utf-8-sig", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
